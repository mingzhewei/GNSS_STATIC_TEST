# GNSS 静态测试分析工具（北云 / 华测）

本仓库包含两套**相互独立但功能一致**的 GNSS 静态测试分析程序，分别针对
**北云（Bynav）** 与 **华测（Huace）** 两款接收机。两套程序在评估算法、
统计口径与单位上保持统一，输出结果具有可比性。

## 程序清单

| 文件 | 对象 | 主数据源 | 说明 |
|------|------|----------|------|
| `gps_rtk_analyzer_bynav.py` | 北云 | ICOM3 中的 **BESTGNSSPOSA**（纯 GNSS）+ **INSPVAXA**（惯导） | 仅分析 COM3；COM4 为二进制流，仅作提示 |
| `gps_rtk_analyzer_huace.py` | 华测 | COM1 中的 **BESTPA**（10 Hz 主数据，BESTPOSA 兜底） | — |
| `gnss_static_hmi.py` | 统一 HMI | 调用上表产品分析脚本 | 产品选择、后台运行、自动打开 HTML 报告 |

## 使用方式

```bash
# 北云
python -X utf8 gps_rtk_analyzer_bynav.py <输入数据文件> -o <输出目录>

# 华测
python -X utf8 gps_rtk_analyzer_huace.py <输入数据文件> -o <输出目录>
```

两套程序均带 GUI（直接运行不带参数即进入界面），并生成 HTML / Markdown 报告。

统一产品选择入口：

```bash
python -X utf8 gnss_static_hmi.py
```

`GNSS Static HMI` 仅作为 HMI/启动器，不包含产品解析逻辑。新增产品时，在
`gnss_static_hmi.py` 的 `PRODUCTS` 注册表中增加一条 `ProductConfig`，指向新的
产品分析脚本；新脚本需兼容 CLI 调用并在输出目录生成 `report.html`。

## 核心准则

1. 北云只使用 **BESTGNSSPOSA**（纯 GNSS）与 **INSPVAXA**（惯导）；
   华测使用 **BESTPA**（10 Hz 主数据，BESTPOSA 仅作兜底）。
2. 弃用 GPIMU；静态测试的精度评估不使用 IMU / 惯导数据，
   INS 相位图 / TRACKSTATA / ENVSTATUSA 等属器件专属章节，予以保留。
3. 全部按**移动站**模式处理；即便输入为运动数据，也按静态测试标准如实评估。
4. 以**实测数据为准**：报文自带时间戳 / 校验位最可信。
5. 对无时间戳、且手册未声明固定周期的报文（如 GSA / GSV），
   周期与均匀性标记为 **N/A（不可实测）**，不伪造统计值。

## 仓库内容说明

本仓库除两套分析程序外，还包含：

- `by_data/`、`huace_data/` —— 实测采集数据（含各自 `*_output/` 分析产物与报告）
- `by_manual/`、`huace_manual/` —— 厂商数据通信接口协议手册（PDF / Markdown / 图片，版权属原厂商，仅作内部参考）

> 说明：上述数据与手册已纳入版本管理。其中两个文件超过 GitHub 建议的 50 MB 阈值（`by_data/_com3.dat` 约 77 MB、`huace_data/...COM11.log` 约 58 MB），推送时远程会提示 LFS 警告，但未超过 100 MB 硬上限，可正常存储。若后续需频繁更新大文件，建议改用 Git LFS 管理。

`.gitignore` 当前仅排除 Python 缓存与系统/临时文件。

## 项目记忆

`.codebuddy/memory/` 存放按日期归档的项目记忆（准则、数据结论、修改记录），
供后续工作延续参考。
