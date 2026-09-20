from __future__ import annotations

import os
import queue
import subprocess
import sys
import threading
import webbrowser
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

APP_TITLE = "GNSS Static HMI"
# 打包成 EXE 后脚本位于 Build 目录；分析器源文件位于其上一级（仓库根目录）。
BASE_DIR = Path(__file__).resolve().parent
ANALYZER_DIR = BASE_DIR.parent if (BASE_DIR.parent / "gps_rtk_analyzer_bynav.py").is_file() else BASE_DIR


def default_output_dir(input_file: str) -> str:
    """当前两个静态分析器的输出目录规则：<输入文件名>_output。"""
    base = os.path.basename(input_file)
    stem = base.rsplit('.', 1)[0] if '.' in base else base
    return os.path.join(os.path.dirname(input_file), stem + '_output')


@dataclass(frozen=True)
class ProductConfig:
    key: str
    display_name: str
    script: str
    filetypes: Sequence[tuple[str, str]]
    color: str = "#34495e"
    extra_args: tuple[str, ...] = ()
    report_name: str = "report.html"
    output_dir_resolver: Callable[[str], str] = default_output_dir


# ============================================================================
# 产品注册表：以后增加新产品，只需要在这里追加一个 ProductConfig
# ============================================================================
PRODUCTS: dict[str, ProductConfig] = {
    "bynav": ProductConfig(
        key="bynav",
        display_name="北云 BYNAV（ICOM3 / COM3 静态）",
        script="gps_rtk_analyzer_bynav.py",
        filetypes=(("北云 COM3 数据", "*.dat"), ("所有文件", "*.*")),
        color="#1a5276",
        extra_args=("--no-open",),
    ),
    "huace": ProductConfig(
        key="huace",
        display_name="华测 HUACE M720（纯 GNSS 静态）",
        script="gps_rtk_analyzer_huace.py",
        filetypes=(("华测日志", "*.dat;*.log"), ("所有文件", "*.*")),
        color="#b03a2e",
        extra_args=("--no-open",),
    ),
}


@dataclass
class AnalysisTask:
    product: ProductConfig
    input_file: str
    status: str = "等待"
    output_dir: str = ""
    report_path: str = ""

    def resolve_output(self) -> None:
        self.output_dir = self.product.output_dir_resolver(self.input_file)
        self.report_path = os.path.join(self.output_dir, self.product.report_name)


def build_command(product: ProductConfig, input_file: str) -> list[str]:
    script_path = ANALYZER_DIR / product.script
    return [sys.executable, "-X", "utf8", str(script_path), input_file, *product.extra_args]


def open_path(path: str) -> None:
    """跨平台打开文件或目录。"""
    if not path:
        return
    if os.name == "nt":
        os.startfile(path)  # type: ignore[attr-defined]
    elif sys.platform == "darwin":
        subprocess.Popen(["open", path])
    else:
        subprocess.Popen(["xdg-open", path])


class GNSSStaticHMI:
    def __init__(self) -> None:
        self.root = tk.Tk()
        self.root.title(APP_TITLE)
        self.root.geometry("1080x720")
        self.root.minsize(920, 620)

        self.tasks: list[AnalysisTask] = []
        self.task_ids: list[int] = []
        self.running = False
        self.worker: threading.Thread | None = None
        self.ui_queue: queue.Queue[tuple] = queue.Queue()

        self._build_ui()
        self.root.after(100, self._poll_ui_queue)

    @classmethod
    def create_for_test(cls) -> "GNSSStaticHMI":
        """测试辅助：不调用 mainloop，仅构造界面对象。"""
        return cls()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def _build_ui(self) -> None:
        self.root.configure(bg="#f5f6fa")

        title = tk.Label(
            self.root,
            text=APP_TITLE,
            font=("微软雅黑", 18, "bold"),
            bg="#f5f6fa",
            fg="#2c3e50",
        )
        title.pack(pady=(18, 6))

        subtitle = tk.Label(
            self.root,
            text="统一调度多个 GNSS 静态产品分析脚本",
            font=("微软雅黑", 10),
            bg="#f5f6fa",
            fg="#7f8c8d",
        )
        subtitle.pack(pady=(0, 14))

        # 产品与文件选择
        select_frame = ttk.LabelFrame(self.root, text="数据源选择")
        select_frame.pack(fill=tk.X, padx=20, pady=(0, 10))

        product_frame = ttk.Frame(select_frame)
        product_frame.pack(fill=tk.X, padx=12, pady=(10, 6))

        ttk.Label(product_frame, text="产品：", width=8).grid(row=0, column=0, sticky=tk.W)
        self.product_var = tk.StringVar(value=next(iter(PRODUCTS.values())).display_name)
        self.product_combo = ttk.Combobox(
            product_frame,
            textvariable=self.product_var,
            state="readonly",
            values=[p.display_name for p in PRODUCTS.values()],
            width=52,
        )
        self.product_combo.grid(row=0, column=1, sticky=tk.W, padx=(0, 12))

        self.browse_button = ttk.Button(product_frame, text="浏览数据文件...", command=self.browse_files)
        self.browse_button.grid(row=0, column=2)

        file_hint = ttk.Label(
            select_frame,
            text="可一次选择多个文件；每个文件会按所选产品加入分析队列。",
            foreground="#7f8c8d",
        )
        file_hint.pack(anchor=tk.W, padx=12, pady=(0, 10))

        # 任务列表
        list_frame = ttk.LabelFrame(self.root, text="分析队列")
        list_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=(0, 10))

        columns = ("product", "file", "status", "output")
        self.task_tree = ttk.Treeview(list_frame, columns=columns, show="headings", selectmode="extended")
        self.task_tree.heading("product", text="产品")
        self.task_tree.heading("file", text="输入文件")
        self.task_tree.heading("status", text="状态")
        self.task_tree.heading("output", text="输出目录")

        self.task_tree.column("product", width=220, anchor=tk.W)
        self.task_tree.column("file", width=300, anchor=tk.W)
        self.task_tree.column("status", width=80, anchor=tk.CENTER)
        self.task_tree.column("output", width=280, anchor=tk.W)

        scrollbar = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=self.task_tree.yview)
        self.task_tree.configure(yscrollcommand=scrollbar.set)
        self.task_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(10, 0), pady=(8, 8))
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y, pady=(8, 8), padx=(0, 10))

        # 控制按钮
        control_frame = ttk.Frame(self.root)
        control_frame.pack(fill=tk.X, padx=20, pady=(0, 10))

        self.start_button = ttk.Button(control_frame, text="开始分析", command=self.start_analysis)
        self.start_button.pack(side=tk.LEFT, padx=(0, 8))

        self.remove_button = ttk.Button(control_frame, text="移除选中", command=self.remove_selected)
        self.remove_button.pack(side=tk.LEFT, padx=(0, 8))

        self.clear_button = ttk.Button(control_frame, text="清空队列", command=self.clear_tasks)
        self.clear_button.pack(side=tk.LEFT, padx=(0, 8))

        self.open_report_button = ttk.Button(control_frame, text="打开选中报告", command=self.open_selected_report)
        self.open_report_button.pack(side=tk.LEFT, padx=(0, 8))

        self.open_dir_button = ttk.Button(control_frame, text="打开选中目录", command=self.open_selected_dir)
        self.open_dir_button.pack(side=tk.LEFT)

        self.open_report_var = tk.BooleanVar(value=True)
        self.open_dir_var = tk.BooleanVar(value=False)

        open_report_check = ttk.Checkbutton(control_frame, text="完成后自动打开 HTML", variable=self.open_report_var)
        open_report_check.pack(side=tk.RIGHT, padx=(8, 0))

        open_dir_check = ttk.Checkbutton(control_frame, text="完成后自动打开目录", variable=self.open_dir_var)
        open_dir_check.pack(side=tk.RIGHT, padx=(8, 0))

        # 日志
        log_frame = ttk.LabelFrame(self.root, text="运行日志")
        log_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=(0, 10))

        self.log_text = tk.Text(log_frame, height=10, wrap=tk.NONE, font=("Consolas", 9))
        log_scroll_y = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        log_scroll_x = ttk.Scrollbar(log_frame, orient=tk.HORIZONTAL, command=self.log_text.xview)
        self.log_text.configure(yscrollcommand=log_scroll_y.set, xscrollcommand=log_scroll_x.set)
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(8, 0), pady=8)
        log_scroll_y.pack(side=tk.RIGHT, fill=tk.Y, pady=8, padx=(0, 8))
        log_scroll_x.pack(side=tk.BOTTOM, fill=tk.X, padx=(8, 8))

        self.log_text.tag_config("info", foreground="#2c3e50")
        self.log_text.tag_config("success", foreground="#27ae60")
        self.log_text.tag_config("error", foreground="#c0392b")

        # 状态栏
        self.status_var = tk.StringVar(value="就绪")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, anchor=tk.W, relief=tk.SUNKEN)
        status_bar.pack(fill=tk.X, side=tk.BOTTOM, padx=20, pady=(0, 8))

    # ------------------------------------------------------------------
    # 队列操作
    # ------------------------------------------------------------------
    def _selected_product(self) -> ProductConfig | None:
        display = self.product_var.get()
        for product in PRODUCTS.values():
            if product.display_name == display:
                return product
        return None

    def browse_files(self) -> None:
        product = self._selected_product()
        if product is None:
            messagebox.showerror(APP_TITLE, "请先选择产品")
            return

        script_path = ANALYZER_DIR / product.script
        if not script_path.is_file():
            messagebox.showerror(APP_TITLE, f"分析脚本不存在：\n{script_path}")
            return

        files = filedialog.askopenfilenames(
            parent=self.root,
            title=f"选择 {product.display_name} 数据文件",
            filetypes=list(product.filetypes),
        )
        if not files:
            return

        added = 0
        for file_path in files:
            if not os.path.isfile(file_path):
                continue
            task = AnalysisTask(product=product, input_file=file_path)
            task.resolve_output()
            self.tasks.append(task)
            self._append_task_to_tree(task)
            added += 1

        if added:
            self._log(f"已添加 {added} 个 {product.display_name} 数据文件。", tag="info")
            self.status_var.set(f"队列中共 {len(self.tasks)} 个任务")
        else:
            messagebox.showwarning(APP_TITLE, "未添加任何有效文件")

    def _append_task_to_tree(self, task: AnalysisTask) -> None:
        item_id = self.task_tree.insert(
            "",
            tk.END,
            values=(task.product.display_name, task.input_file, task.status, task.output_dir),
        )
        self.task_ids.append(item_id)

    def _update_task_status(self, index: int, status: str, report_path: str = "") -> None:
        if index < 0 or index >= len(self.tasks):
            return
        task = self.tasks[index]
        task.status = status
        if report_path:
            task.report_path = report_path
        if index < len(self.task_ids):
            values = (task.product.display_name, task.input_file, task.status, task.output_dir)
            self.task_tree.item(self.task_ids[index], values=values)

    def remove_selected(self) -> None:
        if self.running:
            messagebox.showwarning(APP_TITLE, "分析进行中，不能修改队列")
            return
        selection = self.task_tree.selection()
        if not selection:
            messagebox.showinfo(APP_TITLE, "请先选择要移除的任务")
            return
        indexes = sorted((self.task_ids.index(item) for item in selection), reverse=True)
        for index in indexes:
            self.task_tree.delete(self.task_ids[index])
            del self.task_ids[index]
            del self.tasks[index]
        self.status_var.set(f"队列中共 {len(self.tasks)} 个任务")

    def clear_tasks(self) -> None:
        if self.running:
            messagebox.showwarning(APP_TITLE, "分析进行中，不能清空队列")
            return
        if not self.tasks:
            return
        if not messagebox.askyesno(APP_TITLE, "确定清空当前分析队列？"):
            return
        for item_id in self.task_ids:
            self.task_tree.delete(item_id)
        self.tasks.clear()
        self.task_ids.clear()
        self.status_var.set("队列已清空")

    def open_selected_report(self) -> None:
        selection = self.task_tree.selection()
        if not selection:
            messagebox.showinfo(APP_TITLE, "请先选择一个任务")
            return
        index = self.task_ids.index(selection[0])
        task = self.tasks[index]
        task.resolve_output()
        if not os.path.isfile(task.report_path):
            messagebox.showwarning(APP_TITLE, f"报告不存在：\n{task.report_path}")
            return
        open_path(task.report_path)

    def open_selected_dir(self) -> None:
        selection = self.task_tree.selection()
        if not selection:
            messagebox.showinfo(APP_TITLE, "请先选择一个任务")
            return
        index = self.task_ids.index(selection[0])
        task = self.tasks[index]
        task.resolve_output()
        if not os.path.isdir(task.output_dir):
            messagebox.showwarning(APP_TITLE, f"输出目录不存在：\n{task.output_dir}")
            return
        open_path(task.output_dir)

    # ------------------------------------------------------------------
    # 分析执行
    # ------------------------------------------------------------------
    def start_analysis(self) -> None:
        if self.running:
            messagebox.showwarning(APP_TITLE, "分析正在进行中")
            return
        if not self.tasks:
            messagebox.showwarning(APP_TITLE, "请先添加数据文件")
            return

        if not messagebox.askyesno(APP_TITLE, f"开始分析 {len(self.tasks)} 个任务？"):
            return

        # 在进入后台线程前读取Tk变量，避免跨线程访问Tkinter
        auto_open_report = self.open_report_var.get()
        auto_open_dir = self.open_dir_var.get()

        self.running = True
        self.start_button.configure(state=tk.DISABLED)
        self.browse_button.configure(state=tk.DISABLED)
        self.remove_button.configure(state=tk.DISABLED)
        self.clear_button.configure(state=tk.DISABLED)
        self.status_var.set("分析中...")

        self.worker = threading.Thread(
            target=self._run_tasks,
            kwargs={"auto_open_report": auto_open_report, "auto_open_dir": auto_open_dir},
            daemon=True,
        )
        self.worker.start()

    def _run_tasks(self, auto_open_report: bool = True, auto_open_dir: bool = False) -> None:
        try:
            for index, task in enumerate(self.tasks):
                self.ui_queue.put(("status", index, "运行中", ""))
                self.ui_queue.put(("log", f"\n===== 开始分析: {task.input_file} =====\n", "info"))

                command = build_command(task.product, task.input_file)
                env = os.environ.copy()
                env["PYTHONUNBUFFERED"] = "1"
                env["PYTHONIOENCODING"] = "utf-8"

                creationflags = 0
                if os.name == "nt":
                    creationflags = subprocess.CREATE_NO_WINDOW

                try:
                    process = subprocess.Popen(
                        command,
                        cwd=str(ANALYZER_DIR),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        encoding="utf-8",
                        errors="replace",
                        env=env,
                        creationflags=creationflags,
                        bufsize=1,
                    )
                except Exception as exc:
                    self.ui_queue.put(("status", index, "失败", ""))
                    self.ui_queue.put(("log", f"启动分析进程失败：{exc}\n", "error"))
                    continue

                assert process.stdout is not None
                for line in process.stdout:
                    self.ui_queue.put(("log", line, "info"))
                return_code = process.wait()

                task.resolve_output()
                if return_code == 0 and os.path.isfile(task.report_path):
                    self.ui_queue.put(("status", index, "成功", task.report_path))
                    self.ui_queue.put(("log", f"分析完成：{task.report_path}\n", "success"))
                    if auto_open_report:
                        open_path(task.report_path)
                    if auto_open_dir:
                        open_path(task.output_dir)
                else:
                    self.ui_queue.put(("status", index, "失败", ""))
                    self.ui_queue.put(("log", f"分析失败，返回码 {return_code}\n", "error"))
        finally:
            self.ui_queue.put(("done", auto_open_report))

    def _poll_ui_queue(self) -> None:
        try:
            while True:
                event = self.ui_queue.get_nowait()
                kind = event[0]
                if kind == "log":
                    _, line, tag = event
                    self._log(line, tag=tag)
                elif kind == "status":
                    _, index, status, report_path = event
                    self._update_task_status(index, status, report_path)
                elif kind == "done":
                    _, auto_open_report = event
                    self.running = False
                    self.start_button.configure(state=tk.NORMAL)
                    self.browse_button.configure(state=tk.NORMAL)
                    self.remove_button.configure(state=tk.NORMAL)
                    self.clear_button.configure(state=tk.NORMAL)
                    self.status_var.set("分析完成")
                    success = sum(1 for t in self.tasks if t.status == "成功")
                    failed = sum(1 for t in self.tasks if t.status == "失败")
                    lines = [
                        "分析全部结束。",
                        "",
                        f"任务总数：{len(self.tasks)}",
                        f"成功：{success}",
                        f"失败：{failed}",
                        f"HTML报告：{'已自动打开' if auto_open_report and failed == 0 and success > 0 else '未打开（可在队列中手动打开）'}",
                    ]
                    if failed:
                        messagebox.showwarning(APP_TITLE, "\n".join(lines))
                    else:
                        messagebox.showinfo(APP_TITLE, "\n".join(lines))
        except queue.Empty:
            pass
        self.root.after(100, self._poll_ui_queue)

    def _log(self, message: str, tag: str = "info") -> None:
        self.log_text.insert(tk.END, message, tag)
        self.log_text.see(tk.END)
        # 限制日志长度，避免长时间运行占用过多内存
        line_count = int(self.log_text.index("end-1c").split(".")[0])
        if line_count > 5000:
            self.log_text.delete("1.0", f"{line_count - 5000}.0")

    def run(self) -> None:
        self.root.mainloop()


def main() -> None:
    app = GNSSStaticHMI()
    app.run()


if __name__ == "__main__":
    main()
