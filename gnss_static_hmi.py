#!/usr/bin/env python3
"""GNSS Static HMI: product-selection front end for GNSS static analyzers.

This script intentionally contains no product-specific parsing logic.  It only
selects an input file, launches the registered analyzer as a separate CLI
process, shows its console output, and opens the generated HTML report.

Adding a product
----------------
Add one ``ProductConfig`` entry to ``PRODUCTS`` below.  The analyzer must be a
CLI-compatible Python script that accepts::

    <input_file> -o <output_dir> --no-open

and writes ``report.html`` into the output directory.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import webbrowser


APP_NAME = "GNSS Static HMI"
APP_VERSION = "1.0.0"


@dataclass(frozen=True)
class ProductConfig:
    """Registration record for one product analyzer."""

    key: str
    display_name: str
    script: Path
    filetypes: tuple[tuple[str, str], ...]
    default_extension: str = ".dat"


# ---------------------------------------------------------------------------
# Product registry.  To add a future product, append one ProductConfig entry.
# Keep product parsing and evaluation logic in the product's analyzer script;
# this HMI remains a launcher/report viewer only.
# ---------------------------------------------------------------------------
_BASE_DIR = Path(__file__).resolve().parent

PRODUCTS: dict[str, ProductConfig] = {
    "bynav": ProductConfig(
        key="bynav",
        display_name="北云 Bynav（COM3 / ICOM3）",
        script=_BASE_DIR / "gps_rtk_analyzer_bynav.py",
        filetypes=(
            ("北云 COM3 数据", "*.dat"),
            ("所有文件", "*.*"),
        ),
        default_extension=".dat",
    ),
    "huace": ProductConfig(
        key="huace",
        display_name="华测 Huace（M720 / COM1）",
        script=_BASE_DIR / "gps_rtk_analyzer_huace.py",
        filetypes=(
            ("华测日志", "*.log *.dat"),
            ("所有文件", "*.*"),
        ),
        default_extension=".log",
    ),
}


class AnalyzerProcess:
    """Thread-safe wrapper around a subprocess with non-blocking output reads."""

    def __init__(self, command: list[str], output_dir: Path):
        self.command = command
        self.output_dir = output_dir
        self.process: subprocess.Popen[str] | None = None
        self.output_queue: queue.Queue[str] = queue.Queue()
        self.thread: threading.Thread | None = None

    def start(self) -> None:
        env = os.environ.copy()
        # Ensure readable UTF-8 console output on Windows pipes.
        env["PYTHONIOENCODING"] = "utf-8"
        self.process = subprocess.Popen(
            self.command,
            cwd=str(_BASE_DIR),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            env=env,
        )
        self.thread = threading.Thread(
            target=self._read_output,
            name="analyzer-output-reader",
            daemon=True,
        )
        self.thread.start()

    def _read_output(self) -> None:
        assert self.process is not None and self.process.stdout is not None
        for line in self.process.stdout:
            self.output_queue.put(line.rstrip("\r\n"))
        self.process.stdout.close()
        self.process.wait()
        self.output_queue.put(None)

    def stop(self) -> None:
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()
            try:
                self.process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait(timeout=3)


class GNSSStaticHMI(tk.Tk):
    """Main HMI window."""

    def __init__(self):
        super().__init__()
        self.title(f"{APP_NAME} v{APP_VERSION}")
        self.geometry("920x680")
        self.minsize(820, 600)

        self.product_key = tk.StringVar(value=next(iter(PRODUCTS)))
        self.input_file = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.dpi = tk.StringVar(value="150")

        self.analyzer: AnalyzerProcess | None = None
        self.poll_after_id: str | None = None
        self.last_html_report: Path | None = None

        self._build_style()
        self._build_ui()
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def _build_style(self) -> None:
        style = ttk.Style(self)
        try:
            style.theme_use("vista")
        except tk.TclError:
            pass
        style.configure("Title.TLabel", font=("Microsoft YaHei UI", 16, "bold"))
        style.configure("SubTitle.TLabel", font=("Microsoft YaHei UI", 10))
        style.configure("Analyze.TButton", font=("Microsoft YaHei UI", 10, "bold"))

    def _build_ui(self) -> None:
        outer = ttk.Frame(self, padding=14)
        outer.pack(fill=tk.BOTH, expand=True)

        header = ttk.Frame(outer)
        header.pack(fill=tk.X, pady=(0, 12))
        ttk.Label(header, text=APP_NAME, style="Title.TLabel").pack(anchor=tk.W)
        ttk.Label(
            header,
            text="选择产品与数据源；分析器将在后台运行，完成后自动打开 HTML 报告。",
            style="SubTitle.TLabel",
        ).pack(anchor=tk.W, pady=(2, 0))

        settings = ttk.LabelFrame(outer, text="分析配置", padding=12)
        settings.pack(fill=tk.X)

        ttk.Label(settings, text="产品：").grid(row=0, column=0, sticky=tk.W, pady=6)
        self.product_combo = ttk.Combobox(
            settings,
            textvariable=self.product_key,
            state="readonly",
            values=list(PRODUCTS),
            width=44,
        )
        self.product_combo.grid(row=0, column=1, columnspan=2, sticky=tk.EW, pady=6)
        self.product_combo.bind("<<ComboboxSelected>>", lambda _e: self.update_product_defaults())

        ttk.Button(
            settings,
            text="查看产品配置",
            command=self.show_product_info,
        ).grid(row=0, column=3, sticky=tk.E, padx=(8, 0))

        ttk.Label(settings, text="输入数据文件：").grid(row=1, column=0, sticky=tk.W, pady=6)
        ttk.Entry(settings, textvariable=self.input_file).grid(
            row=1, column=1, columnspan=2, sticky=tk.EW, pady=6
        )
        ttk.Button(settings, text="浏览…", command=self.browse_input).grid(
            row=1, column=3, sticky=tk.E, padx=(8, 0)
        )

        ttk.Label(settings, text="输出目录：").grid(row=2, column=0, sticky=tk.W, pady=6)
        ttk.Entry(settings, textvariable=self.output_dir).grid(
            row=2, column=1, columnspan=2, sticky=tk.EW, pady=6
        )
        ttk.Button(settings, text="选择…", command=self.browse_output).grid(
            row=2, column=3, sticky=tk.E, padx=(8, 0)
        )

        ttk.Label(settings, text="图表 DPI：").grid(row=3, column=0, sticky=tk.W, pady=6)
        self.dpi_spin = ttk.Spinbox(
            settings,
            from_=72,
            to=600,
            increment=10,
            textvariable=self.dpi,
            width=10,
        )
        self.dpi_spin.grid(row=3, column=1, sticky=tk.W, pady=6)

        self.open_output_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            settings,
            text="完成后同时打开输出目录",
            variable=self.open_output_var,
        ).grid(row=3, column=2, sticky=tk.E)

        settings.columnconfigure(1, weight=1)
        settings.columnconfigure(2, weight=1)

        actions = ttk.Frame(outer)
        actions.pack(fill=tk.X, pady=(10, 8))
        self.analyze_button = ttk.Button(
            actions,
            text="分析",
            style="Analyze.TButton",
            command=self.start_analysis,
        )
        self.analyze_button.pack(side=tk.LEFT)
        self.cancel_button = ttk.Button(
            actions,
            text="取消",
            command=self.cancel_analysis,
            state=tk.DISABLED,
        )
        self.cancel_button.pack(side=tk.LEFT, padx=(8, 0))
        self.open_report_button = ttk.Button(
            actions,
            text="打开 HTML 报告",
            command=self.open_last_report,
            state=tk.DISABLED,
        )
        self.open_report_button.pack(side=tk.RIGHT)
        self.open_folder_button = ttk.Button(
            actions,
            text="打开输出目录",
            command=self.open_last_output_dir,
            state=tk.DISABLED,
        )
        self.open_folder_button.pack(side=tk.RIGHT, padx=(0, 8))

        self.status_var = tk.StringVar(value="就绪")
        status_bar = ttk.Label(outer, textvariable=self.status_var, anchor=tk.W)
        status_bar.pack(fill=tk.X, pady=(2, 0))

        log_frame = ttk.LabelFrame(outer, text="分析器输出", padding=8)
        log_frame.pack(fill=tk.BOTH, expand=True, pady=(10, 0))
        self.log_text = tk.Text(
            log_frame,
            wrap=tk.NONE,
            font=("Consolas", 10),
            state=tk.DISABLED,
        )
        y_scroll = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        x_scroll = ttk.Scrollbar(log_frame, orient=tk.HORIZONTAL, command=self.log_text.xview)
        self.log_text.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.log_text.grid(row=0, column=0, sticky=tk.NSEW)
        y_scroll.grid(row=0, column=1, sticky=tk.NS)
        x_scroll.grid(row=1, column=0, sticky=tk.EW)
        log_frame.rowconfigure(0, weight=1)
        log_frame.columnconfigure(0, weight=1)

    # ------------------------------------------------------------------
    # UI helpers
    # ------------------------------------------------------------------
    def selected_product(self) -> ProductConfig:
        return PRODUCTS[self.product_key.get()]

    def update_product_defaults(self) -> None:
        current_input = Path(self.input_file.get()) if self.input_file.get() else None
        if current_input and current_input.is_file():
            self.output_dir.set(str(current_input.with_name(current_input.stem + "_output")))

    def show_product_info(self) -> None:
        p = self.selected_product()
        state = "存在" if p.script.is_file() else "缺失"
        messagebox.showinfo(
            "产品配置",
            "产品注册信息\n\n"
            f"名称：{p.display_name}\n"
            f"注册键：{p.key}\n"
            f"分析脚本：{p.script}\n"
            f"脚本状态：{state}",
        )

    def browse_input(self) -> None:
        p = self.selected_product()
        filename = filedialog.askopenfilename(
            title=f"选择 {p.display_name} 数据文件",
            filetypes=p.filetypes,
        )
        if not filename:
            return
        self.input_file.set(filename)
        input_path = Path(filename)
        self.output_dir.set(str(input_path.with_name(input_path.stem + "_output")))

    def browse_output(self) -> None:
        dirname = filedialog.askdirectory(title="选择输出目录")
        if dirname:
            self.output_dir.set(dirname)

    def append_log(self, line: str) -> None:
        self.log_text.configure(state=tk.NORMAL)
        self.log_text.insert(tk.END, line + "\n")
        self.log_text.see(tk.END)
        self.log_text.configure(state=tk.DISABLED)

    def clear_log(self) -> None:
        self.log_text.configure(state=tk.NORMAL)
        self.log_text.delete("1.0", tk.END)
        self.log_text.configure(state=tk.DISABLED)

    def set_running(self, running: bool) -> None:
        if running:
            self.analyze_button.configure(state=tk.DISABLED)
            self.cancel_button.configure(state=tk.NORMAL)
            self.open_report_button.configure(state=tk.DISABLED)
            self.open_folder_button.configure(state=tk.DISABLED)
        else:
            self.analyze_button.configure(state=tk.NORMAL)
            self.cancel_button.configure(state=tk.DISABLED)
            if self.last_html_report is not None:
                self.open_report_button.configure(state=tk.NORMAL)
            self.open_folder_button.configure(state=tk.NORMAL)

    # ------------------------------------------------------------------
    # Analysis execution
    # ------------------------------------------------------------------
    def validate_inputs(self) -> tuple[Path, Path, int] | None:
        try:
            p = self.selected_product()
        except KeyError:
            messagebox.showerror("错误", "请选择有效产品。")
            return None

        if not p.script.is_file():
            messagebox.showerror("错误", f"产品分析脚本不存在：\n{p.script}")
            return None

        if not self.input_file.get():
            messagebox.showerror("错误", "请选择输入数据文件。")
            return None

        input_path = Path(self.input_file.get())
        if not input_path.is_file():
            messagebox.showerror("错误", f"输入文件不存在：\n{input_path}")
            return None

        output_path = Path(self.output_dir.get()) if self.output_dir.get() else input_path.with_name(
            input_path.stem + "_output"
        )

        try:
            dpi = int(self.dpi.get())
        except ValueError:
            messagebox.showerror("错误", "图表 DPI 必须为整数。")
            return None
        if not 72 <= dpi <= 600:
            messagebox.showerror("错误", "图表 DPI 必须在 72 到 600 之间。")
            return None

        return input_path, output_path, dpi

    def build_command(self, input_path: Path, output_path: Path, dpi: int) -> list[str]:
        p = self.selected_product()
        return [
            sys.executable,
            "-X",
            "utf8",
            str(p.script),
            str(input_path),
            "-o",
            str(output_path),
            "--dpi",
            str(dpi),
            "--no-open",
        ]

    def start_analysis(self) -> None:
        validated = self.validate_inputs()
        if validated is None:
            return
        input_path, output_path, dpi = validated
        command = self.build_command(input_path, output_path, dpi)

        self.last_html_report = None
        self.clear_log()
        self.append_log(f"$ {' '.join(command)}")
        self.status_var.set(f"正在分析：{self.selected_product().display_name}")
        self.set_running(True)

        try:
            output_path.mkdir(parents=True, exist_ok=True)
            self.analyzer = AnalyzerProcess(command, output_path)
            self.analyzer.start()
        except Exception as e:
            self.set_running(False)
            self.status_var.set("启动失败")
            messagebox.showerror("启动失败", f"无法启动分析进程：\n{e}")
            return

        self.poll_after_id = self.after(100, self.poll_analyzer_output)

    def poll_analyzer_output(self) -> None:
        if self.analyzer is None:
            return

        finished = False
        while True:
            try:
                line = self.analyzer.output_queue.get_nowait()
            except queue.Empty:
                break
            if line is None:
                finished = True
                break
            self.append_log(line)

        if finished:
            self.finish_analysis()
            return

        self.poll_after_id = self.after(100, self.poll_analyzer_output)

    def finish_analysis(self) -> None:
        if self.analyzer is None:
            return
        return_code = self.analyzer.process.returncode if self.analyzer.process else -1
        output_dir = self.analyzer.output_dir
        self.analyzer = None
        self.set_running(False)

        if return_code != 0:
            self.status_var.set("分析失败")
            messagebox.showerror(
                "分析失败",
                f"产品分析脚本返回码：{return_code}\n\n请查看下方日志定位问题。",
            )
            return

        html_path = output_dir / "report.html"
        if not html_path.is_file():
            self.status_var.set("分析完成，但未找到报告")
            messagebox.showerror(
                "报告缺失",
                f"分析进程已成功结束，但未找到 HTML 报告：\n{html_path}",
            )
            return

        self.last_html_report = html_path
        self.open_report_button.configure(state=tk.NORMAL)
        self.open_folder_button.configure(state=tk.NORMAL)
        self.status_var.set("分析完成")

        opened = self.open_path_in_browser(html_path)
        if self.open_output_var.get():
            self.open_directory(output_dir)

        state = "已自动打开" if opened else "自动打开失败"
        messagebox.showinfo(
            "分析完成",
            "产品分析已完成。\n\n"
            f"输出目录：{output_dir}\n"
            f"HTML报告：{html_path}\n"
            f"报告状态：{state}",
        )

    def cancel_analysis(self) -> None:
        if self.analyzer is None:
            return
        if not messagebox.askyesno("取消分析", "确定要终止当前分析进程吗？"):
            return
        self.status_var.set("正在取消…")
        self.analyzer.stop()
        # The reader thread will emit None and finish_analysis() will show failure.

    def open_last_report(self) -> None:
        if self.last_html_report is None or not self.last_html_report.is_file():
            messagebox.showerror("错误", "当前没有可打开的 HTML 报告。")
            return
        self.open_path_in_browser(self.last_html_report)

    def open_last_output_dir(self) -> None:
        if self.analyzer is not None:
            self.open_directory(self.analyzer.output_dir)
            return
        # Use the report's parent when no process is running.
        if self.last_html_report is not None:
            self.open_directory(self.last_html_report.parent)
            return
        output = self.output_dir.get()
        if output:
            self.open_directory(Path(output))
            return
        messagebox.showerror("错误", "当前没有可打开的输出目录。")

    @staticmethod
    def open_path_in_browser(path: Path) -> bool:
        try:
            if webbrowser.open(path.resolve().as_uri()):
                return True
            if os.name == "nt":
                os.startfile(str(path))  # type: ignore[attr-defined]
                return True
            if sys.platform == "darwin":
                subprocess.run(["open", str(path)], check=False)
                return True
            subprocess.run(["xdg-open", str(path)], check=False)
            return True
        except Exception:
            return False

    @staticmethod
    def open_directory(path: Path) -> None:
        try:
            if os.name == "nt":
                os.startfile(str(path))  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.run(["open", str(path)], check=False)
            else:
                subprocess.run(["xdg-open", str(path)], check=False)
        except Exception as e:
            messagebox.showerror("错误", f"无法打开目录：\n{e}")

    def on_close(self) -> None:
        if self.analyzer is not None:
            if not messagebox.askokcancel("退出", "分析正在进行中，确定要终止并退出吗？"):
                return
            self.analyzer.stop()
        if self.poll_after_id is not None:
            try:
                self.after_cancel(self.poll_after_id)
            except tk.TclError:
                pass
        self.destroy()


def main() -> int:
    app = GNSSStaticHMI()
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
