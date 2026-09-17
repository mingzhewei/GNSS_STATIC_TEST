# GNSS 静态测试分析工具（北云 / 华测）

本仓库包含两套**相互独立但功能一致**的 GNSS 静态测试分析程序，分别针对
**北云（Bynav）** 与 **华测（Huace）** 两款接收机。两套程序在评估算法、
统计口径与单位上保持统一，输出结果具有可比性。

## 程序清单

| 文件 | 对象 | 主数据源 | 说明 |
|------|------|----------|------|
| `gps_rtk_analyzer_bynav.py` | 北云 | ICOM3 中的 **BESTGNSSPOSA**（纯 GNSS）+ **INSPVAXA**（惯导） | 仅分析 COM3；COM4 为二进制流，仅作提示 |
| `gps_rtk_analyzer_huace.py` | 华测 | COM1 中的 **BESTPA**（10 Hz 主数据，BESTPOSA 兜底） | — |

## 使用方式

```bash
# 北云
python -X utf8 gps_rtk_analyzer_bynav.py <输入数据文件> -o <输出目录>

# 华测
python -X utf8 gps_rtk_analyzer_huace.py <输入数据文件> -o <输出目录>
```

两套程序均带 GUI（直接运行不带参数即进入界面），并生成 HTML / Markdown 报告。

## 核心准则

1. 北云只使用 **BESTGNSSPOSA**（纯 GNSS）与 **INSPVAXA**（惯导）；
   华测使用 **BESTPA**（10 Hz 主数据，BESTPOSA 仅作兜底）。
2. 弃用 GPIMU；静态测试的精度评估不使用 IMU / 惯导数据，
   INS 相位图 / TRACKSTATA / ENVSTATUSA 等属器件专属章节，予以保留。
3. 全部按**移动站**模式处理；即便输入为运动数据，也按静态测试标准如实评估。
4. 以**实测数据为准**：报文自带时间戳 / 校验位最可信。
5. 对无时间戳、且手册未声明固定周期的报文（如 GSA / GSV），
   周期与均匀性标记为 **N/A（不可实测）**，不伪造统计值。

## 不入库的内容

以下内容体积较大或属版权资料，已通过 `.gitignore` 排除，**不纳入版本管理**：

- `by_data/`、`huace_data/` —— 实测采集数据
- `*_output/`、`*_report/` —— 程序生成的分析产物与报告
- `by_manual/`、`huace_manual/`、`*.pdf` —— 厂商数据通信接口协议手册（版权资料）
- `__pycache__/`、`.mypy_cache/` 等缓存

## 项目记忆

`.codebuddy/memory/` 存放按日期归档的项目记忆（准则、数据结论、修改记录），
供后续工作延续参考。
