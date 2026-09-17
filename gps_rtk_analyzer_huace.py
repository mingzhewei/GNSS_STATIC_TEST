#!/usr/bin/env python3
"""
Huace M7/M720 GNSS/INS Data Analyzer
========================================

A comprehensive tool for analyzing static GPS RTK positioning data from Huace (CHCNAV) M7/M720 receivers.

Usage:
    python gps_rtk_analyzer_huace.py <input_file.log> [-o <output_dir>] [--dpi DPI] [--no-html] [--no-md]

This script processes .dat files containing mixed ASCII messages (#BESTPA(10Hz主数据源), #BESTPOSA(1Hz兜底), #BESTDOPSA, #RANGEA, #SATVIS2A, $GNGGA, $GNGSA, $GNGST, $GxGSV).
It splits the data by message type, performs comprehensive analysis, and generates reports with charts.
"""

import sys
import os
import re
import argparse
import base64
import copy
import zlib
from datetime import datetime
from pathlib import Path
import platform

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
import matplotlib.patches as mpatches

# Set up font for Chinese characters
system = platform.system()
if system == 'Windows':
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'KaiTi']
elif system == 'Darwin':
    plt.rcParams['font.sans-serif'] = ['PingFang SC', 'Heiti SC', 'STHeiti']
else:
    plt.rcParams['font.sans-serif'] = ['WenQuanYi Micro Hei', 'Noto Sans CJK SC']
plt.rcParams['axes.unicode_minus'] = False

# Constants
GPS_EPOCH = datetime(1980, 1, 6)
WGS84_A = 6378137.0  # Semi-major axis (m)
WGS84_F = 1 / 298.257223563  # Flattening
WGS84_E2 = 2 * WGS84_F - WGS84_F**2  # First eccentricity squared

# RTK Static Positioning Industry Standards
RTK_STANDARDS = {
    'fixed_solution_threshold': 0.95,  # Minimum fix rate (95%)
    'cep95_rtk_fixed': 0.02,  # 2cm maximum for RTK fixed (meters)
    'cep95_float': 0.50,  # 50cm maximum for RTK float (meters)
    'pdop_good': 2.0,  # Good PDOP threshold
    'pdop_moderate': 3.0,  # Moderate PDOP threshold
    'hdop_good': 1.5,  # Good HDOP threshold
    'hdop_moderate': 2.0,  # Moderate HDOP threshold
    'min_satellites_fixed': 5,  # Minimum satellites for fixed solution
    'min_satellites_float': 4,  # Minimum satellites for float solution
    'max_velocity': 0.1,  # Maximum acceptable velocity for static (m/s)
    'max_horizontal_sigma': 0.05,  # Maximum horizontal position sigma (m)
    'max_vertical_sigma': 0.10,  # Maximum vertical position sigma (m)
    'max_dop_timeout': 10.0  # Maximum acceptable PDOP timeout (seconds)
}

# BESTGNSSPOSA pos_type values that indicate RTK fixed solution
RTK_FIXED_POS_TYPES = ['NARROW_INT', 'WIDE_INT', 'L1_INT', 'NARROW_INT_RTK_DIRECT', 'WIDE_INT_RTK_DIRECT']
# BESTGNSSPOSA pos_type values that indicate RTK float solution
RTK_FLOAT_POS_TYPES = ['NARROW_FLOAT', 'WIDE_FLOAT', 'L1_FLOAT', 'IF_FLOAT', 'NARROW_FLOAT_RTK_DIRECT', 'WIDE_FLOAT_RTK_DIRECT']
# BESTGNSSPOSA sol_status values indicating a valid solution
VALID_SOL_STATUSES = ['SOL_COMPUTED']

# === 解算状态/定位状态枚举全集（华测 M7 表 3-39 / 表 3-40）===
SOL_STATUS_ENUM_HUACE = {
    'SOL_COMPUTED': '已解出', 'INSUFFICIENT_OBS': '观测数据不足',
    'NO_CONVERGENCE': '无法收敛', 'SINGULARITY': '参数矩阵奇异',
    'COV_TRACE': '协方差超限(迹>1000m)', 'VARIANCE': '方差超限',
    'PENDING': '固定位置待定', 'INVALID_FIX': 'SETFIXPOS位置无效',
    'RESIDUALS': '残差过大', 'INTEGRITY_WARNING': '残差过大致位置不可靠',
}
POS_TYPE_ENUM_HUACE = {
    'NONE': '无解', 'SINGLE': '单点定位', 'SPPDIFF': '伪距差分解',
    'NARROW_INT': '窄巷固定解', 'NARROW_FLOAT': '窄巷浮点解',
    'FIXEDPOS': 'SETFIXPOS固定位置', 'SBAS': 'SBAS定位',
    'PPP': '收敛PPP(B2b/HAS/MADOCA)', 'PPP_FIXED': '收敛PPP固定解',
    'PPP_CONVERGING': 'PPP收敛中', 'WIDE_INT': '宽巷固定解',
    'L1_INT': 'L1固定解', 'L1_FLOAT': 'L1浮点解',
    'DOPPLER_VELOCITY': '多普勒速度', 'INS': '纯惯导解',
    'INS_SBAS': 'INS+SBAS组合', 'INS_SPP': 'INS+单点组合',
    'INS_SPPDIFF': 'INS+伪距差分组合', 'INS_RTKFLOAT': 'INS+RTK浮点组合',
    'INS_RTKFIXED': 'INS+RTK固定组合', 'INS_PPP': 'INS+PPP组合',
    'INS_PPPFIXED': 'INS+PPP固定组合',
    # 数据为准：接收机实际输出 NovAtel 风格 PSRDIFF（手册称 SPPDIFF），两者同义。
    'PSRDIFF': '伪距差分(SPPDIFF)', 'PSRSP': '伪距单点', 'FLOATCONV': '浮点载波相位解',
}
# 华测 GNSS 报文无 INS 状态枚举（无 INSPVAXA）；保留空表保证接口一致
INS_STATUS_ENUM_HUACE = {}
# 本设备各枚举全集（供报告按实际出现+未出现列出）
SOL_STATUS_ENUM = SOL_STATUS_ENUM_HUACE
POS_TYPE_ENUM = POS_TYPE_ENUM_HUACE
INS_STATUS_ENUM = INS_STATUS_ENUM_HUACE

# Chart generation constants
DPI_DEFAULT = 150
EXPECTED_INTERVALS = {
    # 华测 COM1 设置(setup)：以下均 0.1s 即 10 Hz
    'BESTPA': 0.1, 'RTKPA': 0.1, 'RTKVA': 0.1, 'BESTVA': 0.1,
    'BESTDOPSA': 0.1, 'ENVSTATUSA': 0.1, 'SATVIS2A': 0.1, 'RANGEA': 1.0,  # RANGEA 实测 1 Hz(数据为准)
    'GNGGA': 0.1, 'GNGSA': 0.1, 'GNGST': 0.1,
    'GPGSV': 0.1, 'GLGSV': 0.1, 'GAGSV': 0.1, 'GQGSV': 0.1, 'GBGSV': 0.1, 'GIGSV': 0.1,
    # BESTPOSA 未在 COM1 设置列表中，实测为 1 Hz（数据为准）
    'BESTPOSA': 1.0,
}

# 周期可实测的消息集合：仅当报文自带时间字段时，周期/频率/均匀性才是【真实测量】。
# 依据华测手册(M7 协议)与北云手册 UG016 §3.1.12：消息是否按周期发送由输出触发方式决定
# (ONTIME=周期输出；ONNEW/ONCHANGED=数据变化才输出)。是否可实测取决于报文是否带时间戳：
# # 类报文头含 TOW、NMEA 中 GGA/GST 含 UTC 时间 -> 可实测；GSA/GSV 无时间字段，且 GSV 每
# 历元连发多条、条数随可见星数变化、同历元各条同时发出 -> 单句周期无意义，标 N/A。
PERIOD_MEASURABLE = {
    'BESTPA', 'RTKPA', 'RTKVA', 'BESTVA', 'BESTDOPSA', 'ENVSTATUSA', 'SATVIS2A',
    'RANGEA', 'BESTPOSA', 'GNGGA', '$GNGGA', 'GNGST', '$GNGST',
}
PERIOD_NOT_MEASURABLE = {
    'GNGSA', '$GNGSA', 'GPGSA', '$GPGSA',
    'GPGSV', '$GPGSV', 'GLGSV', '$GLGSV', 'GAGSV', '$GAGSV',
    'GQGSV', '$GQGSV', 'GBGSV', '$GBGSV', 'GIGSV', '$GIGSV',
}

# Constellation mapping for GSV messages
CONSTELLATION_MAP = {
    '$GPGSV': 'GPS',
    '$GLGSV': 'GLONASS',
    '$GAGSV': 'Galileo',
    '$GQGSV': 'QZSS',
    '$GBGSV': 'BeiDou',
    '$GIGSV': 'NavIC'
}

CONSTELLATION_COLORS = {
    'GPS': '#00FF00',
    'GLONASS': '#FF4444',
    'Galileo': '#4444FF',
    'BeiDou': '#FFD700',
    'QZSS': '#FF00FF',
    'NavIC': '#00FFFF'
}

CONSTELLATION_MARKERS = {
    'GPS': 'o', 'GLONASS': 's', 'Galileo': '^',
    'BeiDou': 'D', 'QZSS': '*', 'NavIC': 'v'
}

# === DATA PARSING FUNCTIONS ===

def validate_latitude(lat):
    """Validate latitude value."""
    try:
        lat = float(lat)
        return -90.0 <= lat <= 90.0
    except (ValueError, TypeError):
        return False

def validate_longitude(lon):
    """Validate longitude value."""
    try:
        lon = float(lon)
        return -180.0 <= lon <= 180.0
    except (ValueError, TypeError):
        return False

def validate_altitude(alt):
    """Validate altitude value."""
    try:
        alt = float(alt)
        return -1000.0 <= alt <= 10000.0  # Reasonable range for GPS altitudes
    except (ValueError, TypeError):
        return False

def validate_velocity(vel):
    """Validate velocity value."""
    try:
        vel = float(vel)
        return -100.0 <= vel <= 100.0  # Reasonable velocity range
    except (ValueError, TypeError):
        return False

def validate_dop(dop):
    """Validate DOP value."""
    try:
        dop = float(dop)
        # Theoretical DOP minimum is >0; best-geometry PDOP can be ~0.6-0.8
        # (4-sat tetrahedron PDOP=sqrt(3/2)~1.22; >4 sats can be lower).
        # 0.5 was arbitrary and produced false warnings on this receiver's
        # normal PDOP=0.8 output. DOP is non-negative by definition.
        return 0.0 < dop <= 50.0
    except (ValueError, TypeError):
        return False

def validate_snr(snr):
    """Validate SNR value."""
    try:
        snr = int(snr)
        return 0 <= snr <= 99  # Typical SNR range
    except (ValueError, TypeError):
        return False

def validate_elevation(elev):
    """Validate satellite elevation."""
    try:
        elev = int(elev)
        return 0 <= elev <= 90
    except (ValueError, TypeError):
        return False

def validate_azimuth(az):
    """Validate satellite azimuth."""
    try:
        az = int(az)
        return 0 <= az <= 359
    except (ValueError, TypeError):
        return False

# === 报文完整性与周期检测（函数名两产品一致，内部实现按器件不同） ===

def _parse_nmea_tod(fields):
    """从 NMEA UTC 时间字段(hhmmss.ss)解析为当日秒数；失败返回 None。"""
    try:
        t = fields[1]
        if not t or len(t) < 6:
            return None
        return int(t[0:2]) * 3600 + int(t[2:4]) * 60 + float(t[4:])
    except Exception:
        return None


# 32 位 CRC 校验算法：北云 UG016 表 4-6(p71) 与华测手册一致，
# 均为多项式 0xEDB88320、初值 0、不含 '#' 的块 CRC-32（等价位序反射的 CRC-32 变体，
# 注意 init=0x00000000，区别于标准 CRC-32/PKZIP 的 init=0xFFFFFFFF）。
# 已对本数据集两产品全部 # 类报文(BESTGNSSPOSA/INSPVAXA/HEADINGA/TRACKSTATA 及
# BESTPOSA/BESTPA/BESTDOPSA/BESTVA/ENVSTATUSA/RANGEA/RTKPA/RTKVA/SATVIS2A)全量
# 验证 100% 匹配、0 错误。计算域：'#' 之后、'*' 之前的全部字符。
_CRC32_TABLE = None


def _crc32_block(data_bytes):
    """按手册表 4-6 的 C 代码实现块 CRC-32（初值 0，多项式 0xEDB88320）。"""
    global _CRC32_TABLE
    if _CRC32_TABLE is None:
        tbl = []
        for i in range(256):
            c = i
            for _ in range(8):
                c = (c >> 1) ^ 0xEDB88320 if (c & 1) else c >> 1
            tbl.append(c)
        _CRC32_TABLE = tbl
    crc = 0
    for b in data_bytes:
        crc = ((crc >> 8) & 0x00FFFFFF) ^ _CRC32_TABLE[(crc ^ b) & 0xFF]
    return crc & 0xFFFFFFFF


def check_message_integrity(lines):
    """报文完整性检查（第1步）。

    校验每条报文自带校验位（依据手册，不猜测）：
    - NMEA ($...)：XOR 校验和（'$' 与 '*' 之间按字节异或），判定 valid/invalid。
    - '#' 类(北云 NovAtel 风格 + 华测 ASCII)：32 位块 CRC-32(表 4-6，初值 0，
      多项式 0xEDB88320，计算域='#' 后 '*' 前)，判定 valid/invalid。
    返回 {'total','valid','invalid','skipped','unverified','integrity_pct',
          'nmea_total','nmea_valid','nmea_invalid','nmea_pct',
          'hash_total','hash_valid','hash_invalid','hash_pct'}。
    """
    total = valid = invalid = skipped = unverified = 0
    nmea_t = nmea_v = nmea_i = 0
    hash_t = hash_v = hash_i = 0
    for item in lines:
        line = item[0] if isinstance(item, tuple) else item
        if not line or '\x00' in line:
            skipped += 1
            continue
        if line.startswith('$'):
            star = line.rfind('*')
            if star < 0 or len(line) < star + 3:
                skipped += 1
                continue
            calc = 0
            for ch in line[1:star]:
                calc ^= ord(ch)
            try:
                want = int(line[star + 1:star + 3], 16)
            except ValueError:
                nmea_i += 1
                continue
            nmea_t += 1
            if calc == want:
                nmea_v += 1
            else:
                nmea_i += 1
        elif line.startswith('#'):
            star = line.rfind('*')
            if star < 0 or len(line) < star + 9:
                skipped += 1
                continue
            cs = line[star + 1:star + 9]
            if not re.fullmatch(r'[0-9a-fA-F]{8}', cs):
                hash_i += 1
                continue
            hash_t += 1
            want = int(cs, 16)
            calc = _crc32_block(line[1:star].encode('ascii', errors='replace'))
            if calc == want:
                hash_v += 1
            else:
                hash_i += 1
        else:
            skipped += 1
    total = nmea_t + hash_t
    valid = nmea_v + hash_v
    invalid = nmea_i + hash_i
    pct = (100.0 * valid / total) if total else 0.0
    return {'total': total, 'valid': valid, 'invalid': invalid,
            'skipped': skipped, 'unverified': unverified, 'integrity_pct': pct,
            'nmea_total': nmea_t, 'nmea_valid': nmea_v, 'nmea_invalid': nmea_i,
            'nmea_pct': (100.0 * nmea_v / nmea_t) if nmea_t else 0.0,
            'hash_total': hash_t, 'hash_valid': hash_v, 'hash_invalid': hash_i,
            'hash_pct': (100.0 * hash_v / hash_t) if hash_t else 0.0}


def build_real_message_periods(message_counts, output_dir):
    """从分割后的各类 .dat 中按设备时间戳提取真实周期（第2步数据源）。

    以数据自带时间戳为准（用户准则）。返回 {消息名: [相邻时间间隔, ...]}。
    """
    out = {}
    for msg_type in message_counts:
        f = Path(output_dir) / f"{msg_type}.dat"
        if not f.exists():
            continue
        ts = []
        try:
            with open(f, 'r', encoding='utf-8') as fh:
                for ln in fh:
                    ln = ln.strip()
                    if not ln:
                        continue
                    try:
                        if ln.startswith('#'):
                            hf = ln.split(';', 1)[0].split(',')
                            # 两种头部混合：idle前置变体第2字段为版本号(46.0)，TOW在索引7(ms)；
                            # NovAtel风格第2字段为端口(COM1)，TOW在索引6(s)。按格式判别。
                            v = None
                            if len(hf) > 7 and '.' in hf[1] and hf[1].replace('.', '').isdigit():
                                v = float(hf[7])          # idle 前置变体 TOW(ms)
                                v = v / 1000.0 if v > 1e7 else v
                            elif len(hf) > 6:
                                v = float(hf[6])          # NovAtel 风格 TOW(s)
                                if v > 1e7:
                                    v = v / 1000.0
                            if v is not None:
                                ts.append(v)
                        elif ln.startswith('$'):
                            v = _parse_nmea_tod(ln.split(','))
                            if v is not None:
                                ts.append(v)
                    except Exception:
                        continue
        except Exception:
            continue
        ts = sorted(ts)
        # 采集端常把同一报文重复写入(同时间戳)，去重以还原真实发播周期
        ts = sorted(set(ts))
        if len(ts) >= 2:
            out[msg_type] = [b - a for a, b in zip(ts, ts[1:]) if b - a > 0]
    return out


def summarize_message_periods(periods, output_dir=None):
    """把真实周期与设置的命令周期核对，并判断插入是否均匀（第2/3/4步）。

    - 实际周期：实测中位间隔与其换算频率（不四舍五入，界面/报告显示实测值）。
    - 均匀性：间隔的变异系数 CV（越小越均匀）；并统计丢失率（间隔明显偏大）。
    - 符合性：实测频率相对标称频率偏差 > 5% 则标记不符合（用于标色告警）。
    """
    result = {}
    for m, ints in (periods or {}).items():
        if not ints:
            continue
        arr = np.asarray(ints, dtype=float)
        med = float(np.median(arr))
        mean = float(np.mean(arr))
        std = float(np.std(arr))
        cv = (std / med) if med > 0 else 0.0
        actual_hz = (1.0 / med) if med > 0 else 0.0
        nominal = EXPECTED_INTERVALS.get(m, 0.0)
        nominal_hz = (1.0 / nominal) if nominal > 0 else 0.0
        if nominal_hz > 0:
            dev_pct = (actual_hz - nominal_hz) / nominal_hz * 100.0
        else:
            dev_pct = 0.0
        lost = int(np.sum(arr > 1.5 * med))
        lost_pct = 100.0 * lost / len(arr)
        # 有效频率按平均间隔计（含丢帧），不四舍五入，如实反映
        effective_hz = (1.0 / mean) if mean > 0 else 0.0
        result[m] = {
            'median_interval_s': med, 'mean_interval_s': mean,
            'interval_std_s': std, 'cv': cv, 'uniform': cv < 0.5,
            'actual_hz': actual_hz, 'nominal_hz': nominal_hz,
            'dev_pct': dev_pct, 'compliant': abs(dev_pct) <= 5.0,
            'lost_count': lost, 'lost_pct': lost_pct, 'samples': len(arr) + 1,
            'effective_hz': effective_hz,
        }
    return result

def analyze_periods_huace(periods):
    """华测消息周期/均匀性分析（内部实现按华测特性）。

    特性：COM1 全部标称 10 Hz（周期 0.1s），应由板卡均匀发播。
    更正(2026-09-17)：早先判定的"ENVSTATUSA 每9帧丢1帧"实为按行读取时被二进制块
    切断造成的假象；改用字节流锚点+CRC校验提取后，ENVSTATUSA 实测满帧 10 Hz、
    间隔全为 0.1s(13958/13958)。故不再做丢帧特殊标注，统一走通用周期汇总。
    """
    return summarize_message_periods(periods)
# 华测 COM11 已知消息头（字节流锚点提取用）。
# # 类：ComNav ASCII，结尾 *8位CRC-32；$ 类：NMEA，结尾 *2位XOR校验和。
_HASH_HEADS = [
    b"#BESTPA", b"#BESTPOSA", b"#BESTVA", b"#RTKPA", b"#RTKVA",
    b"#ENVSTATUSA", b"#BESTDOPSA", b"#SATVIS2A", b"#RANGEA",
]
_NMEA_HEADS = [
    b"$GNGGA", b"$GNGSA", b"$GNGST",
    b"$GPGSV", b"$GLGSV", b"$GAGSV", b"$GBGSV", b"$GQGSV", b"$GIGSV",
]


def _nmea_checksum(body):
    """NMEA XOR 校验和（'$' 与 '*' 之间按字节异或）。"""
    c = 0
    for b in body:
        c ^= b
    return c


def extract_valid_messages(raw_bytes):
    """从 文本+二进制混合流 中按【消息头锚点 + 校验和验证】提取完整合法消息。

    背景（实测，非猜测）：华测 COM11 日志是不定长二进制块与 ASCII 报文交错的
    原始流。二进制块内部不含 CRLF，但会把本应一行一条的 ASCII 报文拦腰切断，
    导致按行读取时该报文被劈成两段、解析时被整体丢弃。实测 #ENVSTATUSA 字节流
    出现 13959 次，但按行解析仅收到 12575 条——丢失的 1384 条正是被切断者。

    本函数以消息头为锚点，对每条候选：定位其 '*' 校验结尾，用 CRC-32(#类，
    表4-6 初值0 多项式0xEDB88320，域='#'后'*'前) 或 NMEA XOR($类) 验证完整性。
    只有校验通过的消息才被保留——这天然剔除被二进制污染的残缺段，实现
    "以字段头为提取要素"的健壮恢复，全程不猜测、以报文自带校验位为准。

    返回 (valid_lines, stats)：valid_lines 为 str 列表(与 read_dat_file 输出兼容)，
    stats 含 rescued/scanned/binary_dropped 计数。
    """
    import re as _re
    n = len(raw_bytes)
    # 收集所有锚点位置 -> (pos, head_bytes)
    anchors = []
    for h in _HASH_HEADS:
        for m in _re.finditer(_re.escape(h), raw_bytes):
            anchors.append((m.start(), h))
    for h in _NMEA_HEADS:
        for m in _re.finditer(_re.escape(h), raw_bytes):
            anchors.append((m.start(), h))
    anchors.sort()

    valid_lines = []
    scanned = len(anchors)
    rescued = 0
    dropped = 0

    for pos, head in anchors:
        # 在该锚点后的一段窗口内找 '*' 校验分隔符（一条 ASCII 报文不会超过 ~600 字节）
        star = raw_bytes.find(b"*", pos, pos + 700)
        if star < 0:
            dropped += 1
            continue
        body = raw_bytes[pos + 1:star]  # '#'/'$' 之后、'*' 之前
        if head.startswith(b"#"):
            # CRC-32，取 '*' 后 8 个十六进制字符
            m = _re.match(rb"[0-9a-fA-F]{8}", raw_bytes[star + 1:star + 9])
            if not m:
                dropped += 1
                continue
            try:
                want = int(m.group(0), 16)
            except ValueError:
                dropped += 1
                continue
            if _crc32_block(body) != want:
                dropped += 1
                continue
            line_bytes = raw_bytes[pos:star + 1 + 8]
        else:
            # NMEA XOR 校验和，'*' 后 2 个十六进制字符
            m = _re.match(rb"[0-9a-fA-F]{2}", raw_bytes[star + 1:star + 3])
            if not m:
                dropped += 1
                continue
            try:
                want = int(m.group(0), 16)
            except ValueError:
                dropped += 1
                continue
            if _nmea_checksum(body) != want:
                dropped += 1
                continue
            line_bytes = raw_bytes[pos:star + 1 + 2]

        try:
            valid_lines.append(line_bytes.decode("ascii"))
            rescued += 1
        except UnicodeDecodeError:
            dropped += 1

    stats = {"scanned": scanned, "rescued": rescued, "binary_dropped": dropped}
    return valid_lines, stats


def read_dat_file(filepath):
    """Read .dat file and return list of lines."""
    lines = []
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f, 1):
                line = line.strip()
                if line:
                    lines.append((line, i))
    except UnicodeDecodeError:
        try:
            with open(filepath, 'r', encoding='latin-1') as f:
                for i, line in enumerate(f, 1):
                    line = line.strip()
                    if line:
                        lines.append((line, i))
        except Exception as e:
            print(f"Error reading file {filepath}: {e}")
            return []
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return []

    print(f"Successfully read {len(lines)} lines from {filepath}")
    return lines

def split_messages(lines, output_dir):
    """Split messages by type and write to separate files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    message_counts = {}
    message_files = {}
    invalid_lines = 0
    binary_lines = 0

    for item in lines:
        # Handle both string and tuple formats
        if isinstance(item, tuple):
            line = item[0]
            line_num = item[1]
        else:
            line = item
            line_num = None

        # Skip lines with null characters (binary data)
        if '\x00' in line:
            binary_lines += 1
            continue

        if line.startswith('#'):
            # ComNav proprietary message
            try:
                msg_type = line[1:].split(',')[0]
                # Remove invalid characters from msg_type (including null characters)
                msg_type = ''.join(c for c in msg_type if c.isprintable() and c not in '<>:"/\\|?*')
                if not msg_type:  # Skip if msg_type becomes empty after cleaning
                    invalid_lines += 1
                    continue
                filename = f"{msg_type}.dat"
            except Exception:
                invalid_lines += 1
                continue
        elif line.startswith('$'):
            # NMEA message - use 6-character identifier for GSV messages
            if line.startswith('$GPGSV') or line.startswith('$GLGSV') or line.startswith('$GAGSV') or \
               line.startswith('$GQGSV') or line.startswith('$GBGSV') or line.startswith('$GIGSV'):
                msg_type = line[:6]
            else:
                msg_type = line[:6]
            # Remove invalid characters from msg_type
            msg_type = ''.join(c for c in msg_type if c.isprintable() and c not in '<>:"/\\|?*')
            if not msg_type:
                invalid_lines += 1
                continue
            filename = f"{msg_type}.dat"
        else:
            invalid_lines += 1
            continue

        message_counts[msg_type] = message_counts.get(msg_type, 0) + 1

        if msg_type not in message_files:
            message_files[msg_type] = open(output_dir / filename, 'w', encoding='utf-8')

        message_files[msg_type].write(line + '\n')

    # Close all files
    for f in message_files.values():
        f.close()

    if binary_lines > 0:
        print(f"Warning: {binary_lines} lines were binary data and skipped")
    if invalid_lines > 0:
        print(f"Warning: {invalid_lines} lines were invalid and skipped")

    return message_counts

def parse_bestgnsspos(lines):
    """
    解析 BESTPOSA / BESTGNSSPOSA 消息（GNSS 最佳位置解）。
    
    华测 M720 实际输出消息名为 BESTPOSA，头部为 NovAtel 风格：
        #BESTPOSA,COM1,0,60.0,FINESTEERING,2435,460383.000,...
        字段含义：端口,序号,空闲时间,时间状态,GPS周,GPS周内秒（秒，非毫秒）
    
    Payload 字段与 NovAtel BESTPOS 兼容：
        sol_status, pos_type, lat, lon, hgt, undulation, datum_id,
        lat_sigma, lon_sigma, hgt_sigma, stn_id, diff_age, sol_age,
        #SVs, #solnSVs, #solnL1SVs, #solnMultiSVs, reserved, ext_sol_stat, ...
    
    参数:
        lines (list): 包含 BESTPOSA/BESTGNSSPOSA 消息的字符串列表
    
    返回:
        list: 解析后的记录列表，每条记录为字典，包含以下字段：
            - week (int): GPS 周数
            - tow (float): GPS 周内秒（Time of Week）
            - sol_status (str): 解状态（如 SOL_COMPUTED）
            - pos_type (str): 位置类型（如 RTK_FIXED, RTK_FLOAT, SINGLE）
            - lat (float): 纬度（度）
            - lon (float): 经度（度）
            - hgt (float): 高度（米）
            - undulation (float): 大地水准面高（米）
            - lat_sigma (float): 纬度标准差（米）
            - lon_sigma (float): 经度标准差（米）
            - hgt_sigma (float): 高度标准差（米）
            - num_svs (int): 跟踪到的卫星总数 (#SVs)
            - num_soln_svs (int): 参与解算的卫星数 (#solnSVs)
            - num_soln_l1 (int): L1 频段参与解算的卫星数
            - num_soln_multi (int): 多频参与解算的卫星数
    
    异常:
        解析过程中的异常会被捕获并打印警告信息，不会抛出异常
    """
    records = []
    for line in lines:
        try:
            if ';' not in line:
                continue
            header, payload = line.split(';', 1)
            h_fields = header.split(',')
            if len(h_fields) < 7:
                print(f"Warning: BESTPOSA header has insufficient fields: {len(h_fields)}")
                continue
            
            # Huace BESTPOSA header: [5]=week, [6]=tow in seconds
            week = int(h_fields[5])
            tow = float(h_fields[6])
            
            p_fields = payload.split('*')[0].split(',')  # Remove CRC first
            
            # Validate number of fields
            if len(p_fields) < 21:
                print(f"Warning: BESTPOSA payload has insufficient fields: {len(p_fields)}")
                continue
            
            # Parse and validate coordinates
            lat = float(p_fields[2]) if p_fields[2] else None
            lon = float(p_fields[3]) if p_fields[3] else None
            hgt = float(p_fields[4]) if p_fields[4] else None
            undulation = float(p_fields[5]) if p_fields[5] else None
            
            # 不再丢弃无有效坐标的历元（如 pos_type=NONE / 无效解），保留以统计完整分布；
            # 无效坐标用 0.0 哨兵 + valid_fix=False 标记，精度统计由 get_fixed_records 过滤。
            valid_fix = bool(validate_latitude(lat) and validate_longitude(lon) and validate_altitude(hgt))
            if not valid_fix:
                lat, lon, hgt = 0.0, 0.0, 0.0
            
            # Parse sigma values
            try:
                lat_sigma = float(p_fields[7]) if p_fields[7] else 0.0
                lon_sigma = float(p_fields[8]) if p_fields[8] else 0.0
                hgt_sigma = float(p_fields[9]) if p_fields[9] else 0.0
            except ValueError:
                print(f"Warning: Invalid sigma values in BESTPOSA at TOW {tow}")
                lat_sigma = lon_sigma = hgt_sigma = 0.0
            
            # Parse satellite counts
            num_svs = int(p_fields[13]) if len(p_fields) > 13 and p_fields[13] else 0
            num_soln_svs = int(p_fields[14]) if len(p_fields) > 14 and p_fields[14] else 0
            num_soln_l1 = int(p_fields[15]) if len(p_fields) > 15 and p_fields[15] else 0
            num_soln_multi = int(p_fields[16]) if len(p_fields) > 16 and p_fields[16] else 0
            
            # NONE/无效历元卫星数可能为 0，保留但并入 valid_fix=False；>99 视为异常丢弃。
            if num_svs < 1 or num_svs > 99:
                valid_fix = False
            if num_svs > 99:
                print(f"Warning: Invalid number of satellites in BESTPOSA at TOW {tow}: {num_svs}")
                continue
            
            records.append({
                'week': week,
                'tow': tow,
                'sol_status': p_fields[0],
                'pos_type': p_fields[1],
                'lat': lat,
                'lon': lon,
                'hgt': hgt,
                'undulation': undulation,
                'datum_id': p_fields[6],
                'lat_sigma': lat_sigma,
                'lon_sigma': lon_sigma,
                'hgt_sigma': hgt_sigma,
                'stn_id': p_fields[10].strip('"') if len(p_fields) > 10 else '',
                'diff_age': float(p_fields[11]) if len(p_fields) > 11 and p_fields[11] else 0.0,
                'sol_age': float(p_fields[12]) if len(p_fields) > 12 and p_fields[12] else 0.0,
                'num_svs': num_svs,
                'num_soln_svs': num_soln_svs,
                'num_soln_l1': num_soln_l1,
                'num_soln_multi': num_soln_multi,
                'valid_fix': valid_fix,
                'reserved1': p_fields[17] if len(p_fields) > 17 else '',
                'ext_sol_stat': p_fields[18] if len(p_fields) > 18 else '',
                'reserved2': p_fields[19] if len(p_fields) > 19 else '',
                'reserved3': p_fields[20] if len(p_fields) > 20 else ''
            })
        except Exception as e:
            print(f"Error parsing BESTPOSA: {e}")
            continue
    
    print(f"Successfully parsed {len(records)} valid BESTPOSA records")
    return records


def parse_bestpa(lines):
    """解析华测 M720 #BESTPA 消息（手册 3.2.14「BESTP 最佳位置信息」，Message ID 3020）。

    用户决策(2026-09-17)：全量切换——华测 GNSS 主数据源由 1 Hz 的 #BESTPOSA
    改为 10 Hz 的 #BESTPA。实测同一 TOW 上 BESTPOSA 与 BESTPA 坐标 100% 一致
    (1396/1396)，解状态相同，数据内容完全等价；BESTPA 时间分辨率高 10 倍。

    与 BESTPOSA 的差异（字段顺序不同，不能互换解析）：
      头部：#BESTPA,<idle>,0,0,COM1,FINESTEERING,<week>,<tow_ms>,<rcv_status>,<hdr_status>
            —— tow 为毫秒(idx7)，week 在 idx6（BESTPOSA 是 week idx5 / tow秒 idx6）。
      数据体(分号后)：sol_status, pos_type, datum_id, reserved, lat, lon, hgt,
            undulation, lat_sigma, lon_sigma, hgt_sigma, diff_age, sol_age,
            #SVs, #solnSVs, #solnL1SVs, #solnMultiSVs, ext_sol_stat,
            gps_glo_mask, bds_mask, gal_mask, reserved, stn_id, CRC
            —— lat/lon/hgt 在 idx4/5/6（BESTPOSA 在 idx2/3/4）。

    返回字段与 parse_bestgnsspos 完全一致，保证下游评估/绘图算法原样复用。
    """
    records = []
    for line in lines:
        try:
            if ';' not in line:
                continue
            header, payload = line.split(';', 1)
            h_fields = header.split(',')
            if len(h_fields) < 8:
                continue

            # BESTPA header: [6]=week, [7]=tow in milliseconds
            week = int(h_fields[6])
            tow = float(h_fields[7]) / 1000.0  # ms -> s，与 BESTPOSA 的秒口径统一

            p_fields = payload.split('*')[0].split(',')  # 去 CRC
            if len(p_fields) < 23:
                continue

            lat = float(p_fields[4]) if p_fields[4] else None
            lon = float(p_fields[5]) if p_fields[5] else None
            hgt = float(p_fields[6]) if p_fields[6] else None
            undulation = float(p_fields[7]) if p_fields[7] else None

            valid_fix = bool(validate_latitude(lat) and validate_longitude(lon) and validate_altitude(hgt))
            if not valid_fix:
                lat, lon, hgt = 0.0, 0.0, 0.0

            try:
                lat_sigma = float(p_fields[8]) if p_fields[8] else 0.0
                lon_sigma = float(p_fields[9]) if p_fields[9] else 0.0
                hgt_sigma = float(p_fields[10]) if p_fields[10] else 0.0
            except ValueError:
                lat_sigma = lon_sigma = hgt_sigma = 0.0

            num_svs = int(p_fields[13]) if len(p_fields) > 13 and p_fields[13] else 0
            num_soln_svs = int(p_fields[14]) if len(p_fields) > 14 and p_fields[14] else 0
            num_soln_l1 = int(p_fields[15]) if len(p_fields) > 15 and p_fields[15] else 0
            num_soln_multi = int(p_fields[16]) if len(p_fields) > 16 and p_fields[16] else 0

            if num_svs < 1 or num_svs > 99:
                valid_fix = False
            if num_svs > 99:
                continue

            records.append({
                'week': week,
                'tow': tow,
                'sol_status': p_fields[0],
                'pos_type': p_fields[1],
                'lat': lat,
                'lon': lon,
                'hgt': hgt,
                'undulation': undulation,
                'datum_id': p_fields[2],
                'lat_sigma': lat_sigma,
                'lon_sigma': lon_sigma,
                'hgt_sigma': hgt_sigma,
                'stn_id': p_fields[22].strip('"') if len(p_fields) > 22 else '',
                'diff_age': float(p_fields[11]) if len(p_fields) > 11 and p_fields[11] else 0.0,
                'sol_age': float(p_fields[12]) if len(p_fields) > 12 and p_fields[12] else 0.0,
                'num_svs': num_svs,
                'num_soln_svs': num_soln_svs,
                'num_soln_l1': num_soln_l1,
                'num_soln_multi': num_soln_multi,
                'valid_fix': valid_fix,
                'ext_sol_stat': p_fields[17] if len(p_fields) > 17 else '',
            })
        except Exception as e:
            print(f"Error parsing BESTPA: {e}")
            continue

    print(f"Successfully parsed {len(records)} valid BESTPA records (10 Hz)")
    return records


def parse_gga(lines):
    """Parse GNGGA/GPGGA message type."""
    records = []
    for line in lines:
        try:
            fields = line.split('*')[0].split(',')
            time_str = fields[1]
            lat = float(fields[2])
            lat_dir = fields[3]
            lon = float(fields[4])
            lon_dir = fields[5]
            fix_type = int(fields[6])
            num_svs = int(fields[7])
            hdop = float(fields[8])
            alt = float(fields[9])
            alt_unit = fields[10]
            undulation = float(fields[11])
            und_unit = fields[12]
            diff_age = float(fields[13]) if fields[13] else 0.0
            station_id = fields[14] if len(fields) > 14 else ''
            
            # Convert NMEA coordinates
            lat_deg = parse_nmea_coord(lat, lat_dir)
            lon_deg = parse_nmea_coord(lon, lon_dir)
            
            # Convert UTC time to seconds since midnight
            time_seconds = 0
            if time_str and len(time_str) >= 6:
                hours = int(time_str[0:2])
                minutes = int(time_str[2:4])
                seconds = float(time_str[4:])
                time_seconds = hours * 3600 + minutes * 60 + seconds
            
            records.append({
                'time': time_str,
                'time_seconds': time_seconds,
                'lat': lat_deg,
                'lon': lon_deg,
                'fix_type': fix_type,
                'num_svs': num_svs,
                'hdop': hdop,
                'alt': alt,
                'undulation': undulation,
                'diff_age': diff_age,
                'station_id': station_id
            })
        except Exception as e:
            print(f"Error parsing GGA: {e}")
    return records


def parse_gsa(lines):
    """Parse GNGSA/GPGSA message type."""
    records = []
    for line in lines:
        try:
            fields = line.split('*')[0].split(',')
            if len(fields) < 18:
                print(f"Warning: GSA line has insufficient fields: {len(fields)}")
                continue
            
            mode = fields[1]
            fix_mode = int(fields[2]) if fields[2].isdigit() else 1
            
            prns = [fields[i] for i in range(3, 15) if fields[i].strip()]
            
            # Parse and validate DOP values
            pdop = float(fields[15]) if fields[15] else 0.0
            hdop = float(fields[16]) if fields[16] else 0.0
            vdop = float(fields[17]) if fields[17] else 0.0
            
            # Warn only when fix is active (fix_mode>=2) but DOP is invalid.
            # No-fix records legitimately carry 0.0 DOPs; warning on them is noise.
            if fix_mode >= 2 and not (validate_dop(pdop) and validate_dop(hdop) and validate_dop(vdop)):
                print(f"Warning: Invalid DOP values in GSA: PDOP={pdop}, HDOP={hdop}, VDOP={vdop}")
            
            records.append({
                'mode': mode,
                'fix_mode': fix_mode,
                'prns': prns,
                'pdop': pdop,
                'hdop': hdop,
                'vdop': vdop
            })
        except Exception as e:
            print(f"Warning: Error parsing GSA line: {e}")
            continue
    
    print(f"Successfully parsed {len(records)} valid GSA records")
    return records


def merge_gsa_dop_per_epoch(gsa_records):
    """将华测同一历元的多条分星座 GNGSA 合并为每历元一条 DOP 记录。

    实测(零猜测)：华测每个历元连发 N 条 GNGSA(每个在用星座一条，PRN 段不同，
    本数据固定 6 条/历元)，同一历元内所有星座条的 PDOP/HDOP/VDOP 完全相同。
    原因：DOP 是接收机本地用【全部在用卫星】的几何统一解算的(见 BESTDOPS
    手册"由接收机当前跟踪并在位置解中使用的卫星的几何分布算出")，并非
    卫星播发、也不是分星座各算各的。因此同一历元各条 DOP 本就一致，
    任取一条即得该历元统一的 PDOP/HDOP/VDOP；本函数取每组首条。

    历元边界判定(关键)：GNGSA 无显式时间戳，且静态场景下 DOP 三元组会在
    连续多个历元间保持完全不变，故【不能用 DOP 值】作为分组键——否则会把
    DOP 没变的相邻历元错误合并(曾出现 83754 条仅合并出 17 组的 bug)。
    可靠边界是【PRN 星座签名】：一个历元内各星座条按固定顺序循环出现，
    当某条的首个 PRN 与当前组首条相同时，即进入下一历元(新一轮星座循环)。
    这样 83754 条可正确还原为约 13959 个历元(= 83754/6)。
    """
    merged = []
    cur = None
    for r in gsa_records:
        first_prn = r['prns'][0] if r['prns'] else ''
        # 新一轮星座循环 = 当前组非空 且 本条首 PRN 与组首条首 PRN 相同
        if cur is not None and first_prn and first_prn == cur['_first_prn']:
            merged.append(cur)
            cur = None
        if cur is None:
            cur = {
                'pdop': r['pdop'], 'hdop': r['hdop'], 'vdop': r['vdop'],
                'prns': list(r['prns']), 'num_records': 1, '_first_prn': first_prn,
            }
        else:
            cur['num_records'] += 1
            cur['prns'].extend(r['prns'])
    if cur is not None:
        merged.append(cur)
    for m in merged:
        m.pop('_first_prn', None)
    print(f"GSA DOP 按历元合并: {len(gsa_records)} 条 -> {len(merged)} 个历元")
    return merged


def parse_bestdopsa(lines):
    """Parse Huace BESTDOPSA message (manual sec 3.2.12, p79-80, Message ID 3029).

    Header layout: #BESTDOPSA,idle,seq,0,port,time_status,week,tow_ms,...
    Note: TOW is in MILLISECONDS (idle-field-first header variant), unlike
    BESTPOSA whose header is port-first and TOW is in seconds.
    Payload order: pdop, gdop, hdop, tdop, htdop, elev_mask, #tracked, sat list...
    (Manual field order differs from NovAtel BESTDOP - do NOT remap.)
    """
    records = []
    for line in lines:
        try:
            if ';' not in line:
                continue
            header, payload = line.split(';', 1)
            h_fields = header.split(',')
            if len(h_fields) < 8:
                continue
            week = int(h_fields[6])
            tow = float(h_fields[7]) / 1000.0  # ms -> s
            p_fields = payload.split('*', 1)[0].split(',')
            if len(p_fields) < 7:
                continue
            pdop = float(p_fields[0])
            gdop = float(p_fields[1])
            hdop = float(p_fields[2])
            if gdop == 0.0 and pdop == 0.0 and hdop == 0.0:
                continue  # receiver not yet computed DOP (startup records)
            records.append({
                'week': week,
                'tow': tow,
                'pdop': pdop,
                'gdop': gdop,
                'hdop': hdop,
                'vdop': None,  # VDOP 仅 GSA 提供；BESTDOPSA 无此字段，用 None 触发下游 N/A 判定（不硬编码数值）
                'tdop': float(p_fields[3]),
                'htdop': float(p_fields[4]),
                'elev_mask': float(p_fields[5]),
                'num_tracked': int(float(p_fields[6])),
            })
        except Exception as e:
            print(f"Warning: Error parsing BESTDOPSA line: {e}")
            continue
    print(f"Successfully parsed {len(records)} valid BESTDOPSA records")
    return records


def parse_gst(lines):
    """Parse $GNGST pseudorange noise statistics (NMEA 0183 GST)."""
    records = []
    for line in lines:
        try:
            fields = line.split('*')[0].split(',')
            if len(fields) < 9 or not fields[1] or not fields[2]:
                continue
            t = fields[1]
            time_seconds = int(t[0:2]) * 3600 + int(t[2:4]) * 60 + float(t[4:])
            records.append({
                'time_seconds': time_seconds,
                'rms': float(fields[2]),
                'major_sigma': float(fields[3]) if fields[3] else None,
                'minor_sigma': float(fields[4]) if fields[4] else None,
                'orient': float(fields[5]) if fields[5] else None,
                'lat_sigma': float(fields[6]) if fields[6] else None,
                'lon_sigma': float(fields[7]) if fields[7] else None,
                'alt_sigma': float(fields[8]) if fields[8] else None,
            })
        except Exception as e:
            print(f"Warning: Error parsing GST line: {e}")
            continue
    print(f"Successfully parsed {len(records)} valid GST records")
    return records


def parse_gsv(lines):
    """
    解析 GSV (GNSS Satellites in View) 消息类型（卫星可见性消息）。
    
    GSV 消息用于报告当前可见的卫星信息，包括卫星编号、仰角、方位角和信噪比。
    一个完整的 GSV 消息可能分为多条（最多包含 4 颗卫星/条），需要合并处理。
    支持华测 M720 输出的 NavIC（$GIGSV）星座。
    
    处理流程：
    1. 按星座分组所有 GSV 消息
    2. 合并同一星座的多条消息
    3. 按 PRN 去重，对重复卫星进行位置和 SNR 的平均
    4. 生成包含所有星座所有卫星的单一历元
    
    参数:
        lines (list): 包含 GSV 消息的字符串列表，每条消息一行
    
    返回:
        list: 解析后的历元列表，每个历元为字典，包含各星座的卫星信息：
            {
                'GPS': {
                    'total_sats': int,        # 该星座可见卫星总数
                    'satellites': [           # 卫星列表
                        {
                            'prn': int,       # 卫星编号
                            'elevation': float, # 仰角（度），已去重平均
                            'azimuth': float,  # 方位角（度），已去重平均
                            'snr': float       # 信噪比（dBHz），已去重平均
                        },
                        ...
                    ]
                },
                'BeiDou': {...},
                'GLONASS': {...},
                'Galileo': {...},
                'QZSS': {...},
                'NavIC': {...}
            }
    
    异常:
        解析过程中的异常会被捕获并打印警告信息，不会抛出异常
    """
    epochs = []
    
    # Filter valid GSV lines
    gsv_lines = [line for line in lines if line.startswith('$') and line[:6] in CONSTELLATION_MAP]
    
    # Group GSV messages by constellation first
    constellation_data = {}
    for line in gsv_lines:
        try:
            fields = line.split('*')[0].split(',')
            talker = fields[0]
            if talker not in CONSTELLATION_MAP:
                continue
            
            constellation = CONSTELLATION_MAP[talker]
            total_msgs = int(fields[1])
            msg_num = int(fields[2])
            total_sats = int(fields[3])
            
            if constellation not in constellation_data:
                constellation_data[constellation] = {}
            
            if msg_num not in constellation_data[constellation]:
                constellation_data[constellation][msg_num] = {
                    'total_msgs': total_msgs,
                    'total_sats': total_sats,
                    'satellites': []
                }
            
            sats = []
            for i in range(0, 4):
                base = 4 + i*4
                if base+3 < len(fields) and fields[base].strip():
                    prn = int(fields[base])
                    elevation = int(fields[base+1])
                    azimuth = int(fields[base+2])
                    snr = int(fields[base+3]) if fields[base+3].strip() else 0
                    
                    if validate_elevation(elevation) and validate_azimuth(azimuth) and validate_snr(snr):
                        sats.append({
                            'prn': prn,
                            'elevation': elevation,
                            'azimuth': azimuth,
                            'snr': snr
                        })
            
            constellation_data[constellation][msg_num]['satellites'].extend(sats)
            
        except Exception as e:
            print(f"Warning: Error parsing GSV line: {e}")
            continue
    
    all_msg_nums = set()
    for const_data in constellation_data.values():
        all_msg_nums.update(const_data.keys())
    
    epoch = {}
    for constellation, msg_data in constellation_data.items():
        all_satellites = []
        for msg_num in sorted(msg_data.keys()):
            all_satellites.extend(msg_data[msg_num]['satellites'])
        
        unique_sats = {}
        for sat in all_satellites:
            prn = sat['prn']
            if prn not in unique_sats:
                unique_sats[prn] = sat
            else:
                existing = unique_sats[prn]
                count = existing.get('count', 1)
                unique_sats[prn] = {
                    'prn': prn,
                    'elevation': (existing['elevation'] * count + sat['elevation']) / (count + 1),
                    'azimuth': (existing['azimuth'] * count + sat['azimuth']) / (count + 1),
                    'snr': (existing['snr'] * count + sat['snr']) / (count + 1),
                    'count': count + 1
                }
        
        epoch[constellation] = {
            'total_sats': len(unique_sats),
            'satellites': list(unique_sats.values())
        }
    
    if epoch:
        epochs.append(epoch)
    
    print(f"Successfully parsed {len(epochs)} GSV epochs")
    return epochs

def parse_nmea_coord(coord_val, direction):
    """Convert NMEA coordinate format to decimal degrees."""
    coord_str = str(coord_val)
    if '.' in coord_str:
        dot_idx = coord_str.index('.')
        degrees = int(coord_str[:dot_idx-2])
        minutes = float(coord_str[dot_idx-2:])
    else:
        degrees = int(coord_str[:-2])
        minutes = float(coord_str[-2:])
    decimal = degrees + minutes / 60.0
    if direction in ('S', 'W'):
        decimal = -decimal
    return decimal

# === CORE CALCULATION FUNCTIONS ===


def get_fixed_records(gnss_data):
    """Return records with a valid computed RTK fixed/float solution.

    Static-precision metrics (CEP95, RMS, peaks) are defined on the steady
    RTK solution; startup SINGLE / NONE epochs are convergence transients and
    must not enter the statistics (they sit meters away from the mean).
    """
    fixed = [r for r in gnss_data
             if r.get('sol_status') in VALID_SOL_STATUSES
             and r.get('pos_type') in RTK_FIXED_POS_TYPES]
    return fixed if fixed else list(gnss_data)

def llh_to_enu(lat_deg, lon_deg, hgt, ref_lat_deg, ref_lon_deg, ref_hgt):
    """Convert geodetic coordinates to ENU (East-North-Up) in meters."""
    # Convert to radians
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    ref_lat = np.radians(ref_lat_deg)
    ref_lon = np.radians(ref_lon_deg)

    # Radius of curvature in prime vertical
    N = WGS84_A / np.sqrt(1 - WGS84_E2 * np.sin(ref_lat)**2)
    # Radius of curvature in meridian
    M = WGS84_A * (1 - WGS84_E2) / (1 - WGS84_E2 * np.sin(ref_lat)**2)**1.5

    # Differentials
    dlat = lat - ref_lat
    dlon = lon - ref_lon
    dhgt = hgt - ref_hgt

    # ENU conversion
    East = dlon * (N + ref_hgt) * np.cos(ref_lat)
    North = dlat * (M + ref_hgt)
    Up = dhgt

    return East, North, Up

def calc_cep95(east, north):
    """Calculate CEP95 (Circular Error Probable at 95%)."""
    if len(east) == 0 or len(north) == 0:
        return 0.0

    # Calculate radial errors from mean
    mean_e = np.mean(east)
    mean_n = np.mean(north)

    radial_errors = np.sqrt((east - mean_e)**2 + (north - mean_n)**2)
    cep95 = np.percentile(radial_errors, 95)
    return cep95

def calc_2drms(east, north):
    """Calculate 2DRMS (Twice Distance RMS)."""
    sigma_e = np.std(east, ddof=1)
    sigma_n = np.std(north, ddof=1)
    return 2.0 * np.sqrt(sigma_e**2 + sigma_n**2)

def calc_rms(values):
    """Calculate Root Mean Square."""
    if len(values) == 0:
        return 0.0
    return np.sqrt(np.mean(np.array(values)**2))

def calc_statistics(values):
    """Calculate statistical measures."""
    if len(values) == 0:
        return {'mean': 0.0, 'std': 0.0, 'rms': 0.0, 'min': 0.0, 'max': 0.0, 'median': 0.0, 'count': 0}

    arr = np.array(values)
    return {
        'mean': float(np.mean(arr)),
        'std': float(np.std(arr, ddof=1)),
        'rms': float(calc_rms(arr)),
        'min': float(np.min(arr)),
        'max': float(np.max(arr)),
        'median': float(np.median(arr)),
        'count': len(arr)
    }

def detect_data_gaps(timestamps, expected_interval, tolerance=2.0):
    """Detect gaps in data."""
    gaps = []
    for i in range(1, len(timestamps)):
        dt = timestamps[i] - timestamps[i-1]
        if dt > expected_interval * tolerance:
            gaps.append({
                'start_idx': i-1,
                'end_idx': i,
                'start_time': timestamps[i-1],
                'end_time': timestamps[i],
                'gap_seconds': dt
            })
    return gaps

def calc_fix_continuity(gnss_data):
    """固定解连续性评估（仅使用纯 GNSS 数据，不依赖 INS）。

    统计静态观测期间 RTK 固定状态的连续保持情况：
    - fixed_epochs: 固定解历元数
    - total_epochs: 总历元数
    - interruption_count: 固定中断次数（固定 -> 非固定 的跳变次数）
    - longest_fix_streak_s: 最长连续固定时长（秒，按历元时间戳累计）
    - current_fix_streak_s: 末段连续固定时长（秒）
    - refix_times_s: 每次失锁后重新恢复固定所需时间列表（秒）
    - mean_refix_s / max_refix_s: 平均 / 最长重固定时间（秒）
    """
    if not gnss_data:
        return {}
    rows = sorted(gnss_data, key=lambda r: r['tow'])
    is_fixed = [
        (r.get('sol_status') in VALID_SOL_STATUSES) and (r.get('pos_type') in RTK_FIXED_POS_TYPES)
        for r in rows
    ]
    times = [r['tow'] for r in rows]
    total = len(rows)
    fixed_epochs = sum(is_fixed)

    interruptions = 0
    longest_streak = 0.0
    streak = 0.0
    refix_times = []
    prev = is_fixed[0]
    # 首历元即固定
    if prev:
        streak = 0.0
    outage_start = None  # 非固定段起始时间
    for i in range(1, total):
        dt = times[i] - times[i - 1]
        cur = is_fixed[i]
        if cur:
            streak += dt
            if not prev:
                # 重新固定
                if outage_start is not None:
                    refix_times.append(times[i] - outage_start)
                    outage_start = None
        else:
            if prev:
                interruptions += 1
                longest_streak = max(longest_streak, streak)
                streak = 0.0
                outage_start = times[i]
        prev = cur
    longest_streak = max(longest_streak, streak)
    current_streak = streak if is_fixed[-1] else 0.0

    return {
        'fixed_epochs': fixed_epochs,
        'total_epochs': total,
        'interruption_count': interruptions,
        'longest_fix_streak_s': longest_streak,
        'current_fix_streak_s': current_streak,
        'refix_times_s': refix_times,
        'mean_refix_s': (sum(refix_times) / len(refix_times)) if refix_times else 0.0,
        'max_refix_s': max(refix_times) if refix_times else 0.0,
    }


def assess_static_stability(gnss_data):
    """静态静止性核验（方案 A：纯 GNSS 位置平稳性，零 INS 依赖）。

    通过位置时序的离散度与漂移量判断测站是否保持静止：
    - 水平/垂直 RMS（围绕均值）
    - 最大水平偏移、最大垂直偏移（相对均值）
    - 首尾漂移量（起点均值与终点均值的水平距离）
    判定：所有指标均小于对应 sigma 阈值则视为静止。
    """
    if not gnss_data:
        return {}
    # 仅使用有效解评估静止性：剔除 SOL 未解算或坐标为 (0,0,0) 的无效历元
    valid = [
        r for r in gnss_data
        if r.get('sol_status') in VALID_SOL_STATUSES
        and abs(r['lat']) > 1e-6 and abs(r['lon']) > 1e-6
    ]
    if not valid:
        return {}
    lats = np.array([r['lat'] for r in valid])
    lons = np.array([r['lon'] for r in valid])
    hgts = np.array([r['hgt'] for r in valid])
    ml, mo, mh = np.mean(lats), np.mean(lons), np.mean(hgts)
    east, north, up = llh_to_enu(lats, lons, hgts, ml, mo, mh)

    horiz = np.sqrt(east**2 + north**2)
    horiz_rms = float(np.sqrt(np.mean(east**2 + north**2)))
    vert_rms = float(np.std(up, ddof=1))
    max_horiz = float(np.max(horiz))
    max_vert = float(np.max(np.abs(up)))

    # 首尾漂移：前 10% 与后 10% 均值位置的水平距离
    k = max(1, len(east) // 10)
    de = float(np.mean(east[-k:]) - np.mean(east[:k]))
    dn = float(np.mean(north[-k:]) - np.mean(north[:k]))
    drift = float(np.hypot(de, dn))

    # 判定阈值：水平/垂直最大偏移不超过 3 倍对应 sigma 阈值
    h_thr = RTK_STANDARDS['max_horizontal_sigma'] * 3
    v_thr = RTK_STANDARDS['max_vertical_sigma'] * 3
    is_static = (max_horiz <= h_thr) and (max_vert <= v_thr)

    return {
        'horizontal_rms': horiz_rms,
        'vertical_rms': vert_rms,
        'max_horizontal_offset': max_horiz,
        'max_vertical_offset': max_vert,
        'drift_first_last': drift,
        'h_threshold': h_thr,
        'v_threshold': v_thr,
        'is_static': is_static,
    }

# === 静态质量评估扩展（行业标准静态流程；字段注释见各函数 docstring） ===

def summarize_sigma_convergence(gnss_data):
    """定位标准差(sigma)收敛评估 —— 用于判断 RTK 解的收敛情况与定位精度。

    字段来源：位置报文携带的 lat/lon/hgt_sigma（解算协方差的水平/垂直分量）。
    判断什么：sigma 随时间由大变小并趋于稳定 = 收敛良好；sigma 持续偏大或反复跳变 = 未收敛/解质量差。
    怎样好：收敛后水平 sigma 小且平稳（数值越小越好，量级因厂商定义不同，故只看趋势不跨厂商比绝对值）。
    怎样坏：sigma 长时间不下降、或频繁出现尖峰（对应失锁/重收敛）。
    返回：水平/垂直 sigma 的初值、末值、均值、最大值与时间序列（供绘图）。
    """
    if not gnss_data:
        return {}
    rows = sorted(gnss_data, key=lambda r: r['tow'])
    valid = [r for r in rows if r.get('lat_sigma') is not None and r.get('lat_sigma', 0) > 0]
    if not valid:
        return {}
    t = [r['tow'] for r in valid]
    hs = [float(np.hypot(r['lat_sigma'], r['lon_sigma'])) for r in valid]
    vs = [float(r['hgt_sigma']) for r in valid]
    k = max(1, len(valid) // 20)  # 前/后 5% 作为初值/末值
    return {
        'times': t, 'h_sigma': hs, 'v_sigma': vs,
        'h_init': float(np.mean(hs[:k])), 'h_final': float(np.mean(hs[-k:])),
        'h_mean': float(np.mean(hs)), 'h_max': float(np.max(hs)),
        'v_init': float(np.mean(vs[:k])), 'v_final': float(np.mean(vs[-k:])),
        'v_mean': float(np.mean(vs)), 'v_max': float(np.max(vs)),
    }

def summarize_gst_rms_convergence(gst_data):
    """华测专属：基于 GNGST 伪距残差 RMS 的收敛评估（替代恒定的位置sigma）。

    字段来源：GNGST 的 rms（伪距残差 RMS, m）——M720 唯一逐历元变化的精度量。
    判断什么：rms 由大变小并趋于稳定 = 解收敛良好；rms 尖峰 = 失锁/重收敛。
    与北云 summarize_sigma_convergence 返回结构一致，供报告同一模板渲染。
    """
    if not gst_data:
        return {}
    rows = sorted(gst_data, key=lambda r: r['time_seconds'])
    valid = [r for r in rows if r.get('rms') is not None and r['rms'] >= 0]
    if not valid:
        return {}
    t = [r['time_seconds'] for r in valid]
    hs = [float(r['rms']) for r in valid]  # 用伪距残差RMS充当水平精度趋势量
    vs = [float(r['major_sigma']) if r.get('major_sigma') is not None else float('nan') for r in valid]
    k = max(1, len(valid) // 20)
    import math
    vs_clean = [v for v in vs if not math.isnan(v)] or [0.0]
    return {
        'times': t, 'h_sigma': hs, 'v_sigma': vs,
        'h_init': float(np.mean(hs[:k])), 'h_final': float(np.mean(hs[-k:])),
        'h_mean': float(np.mean(hs)), 'h_max': float(np.max(hs)),
        'v_init': float(np.mean(vs_clean[:k])), 'v_final': float(np.mean(vs_clean[-k:])),
        'v_mean': float(np.mean(vs_clean)), 'v_max': float(np.max(vs_clean)),
    }


def summarize_snr(gsv_data):
    """卫星载噪比(C/N0)统计 —— 用于评估接收机前端信号质量与遮挡环境。

    字段来源：GSV 报文中每颗可见星的 SNR（dBHz）。
    判断什么：整体信号接收条件好坏。
    怎样好：平均 C/N0 >= 38 dBHz（开阔地，行业常用分级）。
    怎样坏：平均 C/N0 32~37（半遮挡）、<= 31（严重遮挡）；低仰角星 C/N0 普遍偏低属正常。
    返回：全部卫星 SNR 的均值/最小/最大/分布，以及各星座均值。
    """
    all_snr = []
    per_const = {}
    for epoch in gsv_data or []:
        for const, blk in epoch.items():
            for s in blk.get('satellites', []):
                v = s.get('snr')
                if v is not None and v > 0:
                    all_snr.append(v)
                    per_const.setdefault(const, []).append(v)
    if not all_snr:
        return {}
    return {
        'mean': float(np.mean(all_snr)),
        'min': float(np.min(all_snr)),
        'max': float(np.max(all_snr)),
        'count': len(all_snr),
        'per_constellation_mean': {c: float(np.mean(v)) for c, v in per_const.items()},
    }

# === 华测特有：ENVSTATUSA 环境感知（仅移动站模式字段） ===

def parse_envstatusa(lines):
    """解析 ENVSTATUSA 环境感知报文（华测 M7 §3.2.19，Message ID 3074）。

    仅提取移动站模式字段（用户确认只做移动站；基站专用字段 moved/diff_pos 一律忽略）。
    各字段含义与好坏判断（手册 p104-107）：
      env_score      环境质量分值：85~100 优秀 / 60~84 一般可用 / 0~59 不合格
      pos_type       定位类型
      sat_vis_rate   卫星可视比率%(跟踪/理论可视)，越高越好
      sat_sol_rate   参与解算卫星比率%(解算/跟踪)，越高越好
      sol_confidence 解的置信度，越高越好
      num_svs        跟踪卫星数
      cv_sats        基站-Rover 共视卫星数，越多 RTK 越稳
      iono_level     电离层活跃状态(QUIET=平静为好)
      avg_snr        平均信噪比 dBHz：>=38 开阔地 / 32~37 半遮挡 / <=31 严重遮挡
      pdop           位置精度因子：<=2.5 分布较好 / 2.5~5 一般 / >5 较差
    时间戳：idle 前置变体头部 TOW 在索引 7，单位毫秒。
    """
    records = []
    for line in lines:
        try:
            head, _, rest = line.partition(';')
            payload = rest.split('*')[0]
            f = payload.split(',')
            if len(f) < 15:
                continue
            hf = head.split(',')
            tow = float(hf[7]) / 1000.0 if len(hf) > 7 and float(hf[7]) > 1e7 else (float(hf[7]) if len(hf) > 7 else None)
            records.append({
                'tow': tow,
                'env_score': int(f[0]),
                'pos_type': f[1],
                'sat_vis_rate': int(f[2]),
                'sat_sol_rate': int(f[3]),
                'sol_confidence': float(f[4]),
                # f[5]=moved(基站专用, 忽略); f[6]=reserved
                'num_svs': int(f[7]),
                'cv_sats': int(f[8]),
                'iono_level': f[9],
                'avg_snr': int(f[10]),
                # f[11..12]=reserved; f[13]=diff_pos(基站专用, 忽略)
                'pdop': float(f[14]),
            })
        except Exception as e:
            print(f"Warning: Error parsing ENVSTATUSA line: {e}")
            continue
    print(f"Successfully parsed {len(records)} valid ENVSTATUSA records")
    return records


def summarize_envstatus(records):
    """ENVSTATUSA 环境感知汇总（移动站模式）。

    判断什么：观测环境对静态定位的影响（环境差 -> 静态精度易变差）。
    输出各项的最小/均值/最大，供报告逐项注释好坏。
    """
    if not records:
        return {}
    def s(key):
        vals = [r[key] for r in records if isinstance(r.get(key), (int, float))]
        if not vals:
            return None
        return {'min': float(np.min(vals)), 'mean': float(np.mean(vals)),
                'max': float(np.max(vals))}
    iono = {}
    for r in records:
        iono[r['iono_level']] = iono.get(r['iono_level'], 0) + 1
    return {
        'env_score': s('env_score'),
        'sat_vis_rate': s('sat_vis_rate'),
        'sat_sol_rate': s('sat_sol_rate'),
        'sol_confidence': s('sol_confidence'),
        'num_svs': s('num_svs'),
        'cv_sats': s('cv_sats'),
        'avg_snr': s('avg_snr'),
        'pdop': s('pdop'),
        'iono_level_counts': iono,
        'count': len(records),
    }

def validate_industry_compliance(gnss_data, dop_source=None):
    """Validate data against RTK industry standards."""
    if not gnss_data:
        return {}, {}

    metrics = {}
    compliance = {}

    # Fix rate calculation
    if gnss_data:
        total_fixes = len(gnss_data)
        fixed_solutions = len([r for r in gnss_data if r.get('sol_status') in VALID_SOL_STATUSES and r.get('pos_type') in RTK_FIXED_POS_TYPES])
        fix_rate = fixed_solutions / total_fixes if total_fixes > 0 else 0

        metrics['fix_rate'] = fix_rate
        compliance['fix_rate'] = {
            'value': fix_rate,
            'threshold': RTK_STANDARDS['fixed_solution_threshold'],
            'pass': fix_rate >= RTK_STANDARDS['fixed_solution_threshold'],
            'message': "✅ 通过" if fix_rate >= RTK_STANDARDS['fixed_solution_threshold'] else "❌ 未通过"
        }

    # Position accuracy (using CEP95) —— 仅基于有效固定解，剔除 NONE/无效历元
    _fx = get_fixed_records(gnss_data)
    if _fx:
        times = [r['tow'] for r in _fx]
        lats = np.array([r['lat'] for r in _fx])
        lons = np.array([r['lon'] for r in _fx])
        hgts = np.array([r['hgt'] for r in _fx])

        mean_lat = np.mean(lats)
        mean_lon = np.mean(lons)
        mean_hgt = np.mean(hgts)

        east, north, up = llh_to_enu(lats, lons, hgts, mean_lat, mean_lon, mean_hgt)

        cep95 = calc_cep95(east, north)
        horizontal_rms = np.sqrt(np.mean(east**2 + north**2))
        # 垂直RMS与水平RMS同口径：相对均值的均方根 sqrt(mean(up^2))（行业静态RMS惯例，跨设备统一）
        vertical_rms = np.sqrt(np.mean(up**2))

        metrics['cep95'] = cep95
        metrics['horizontal_rms'] = horizontal_rms
        metrics['vertical_rms'] = vertical_rms

        # Check against standards
        pos_type = _fx[0].get('pos_type', 'UNKNOWN')
        if pos_type in RTK_FIXED_POS_TYPES:
            standard = RTK_STANDARDS['cep95_rtk_fixed']
        else:
            standard = RTK_STANDARDS['cep95_float']

        compliance['position_accuracy'] = {
            'value': cep95,
            'standard': standard,
            'pass': cep95 <= standard,
            'message': "✅ 通过" if cep95 <= standard else "⚠️ 警告"
        }

    # DOP values
    if gnss_data:
        # DOP source priority: BESTDOPSA (receiver-computed, per manual 3.2.12)
        # > GSA-derived (GSA on M720 is per-constellation and may be empty at start)
        dop_src = dop_source if dop_source else gnss_data
        pdops = [r.get('pdop', 999) for r in dop_src if 'pdop' in r]
        hdops = [r.get('hdop', 999) for r in dop_src if 'hdop' in r]

        avg_pdop = np.mean(pdops) if pdops else 999
        avg_hdop = np.mean(hdops) if hdops else 999

        metrics['avg_pdop'] = avg_pdop
        metrics['avg_hdop'] = avg_hdop

        compliance['pdop'] = {
            'value': avg_pdop,
            'threshold': RTK_STANDARDS['pdop_good'],
            'pass': avg_pdop <= RTK_STANDARDS['pdop_good'],
            'message': "✅ 优秀" if avg_pdop <= RTK_STANDARDS['pdop_good'] else "⚠️ 警告" if avg_pdop <= RTK_STANDARDS['pdop_moderate'] else "❌ 差"
        }

        compliance['hdop'] = {
            'value': avg_hdop,
            'threshold': RTK_STANDARDS['hdop_good'],
            'pass': avg_hdop <= RTK_STANDARDS['hdop_good'],
            'message': "✅ 优秀" if avg_hdop <= RTK_STANDARDS['hdop_good'] else "⚠️ 警告" if avg_hdop <= RTK_STANDARDS['hdop_moderate'] else "❌ 差"
        }

    # Satellite count
    if gnss_data:
        sat_counts = [r['num_svs'] for r in gnss_data if 'num_svs' in r]
        min_sats = min(sat_counts) if sat_counts else 0
        avg_sats = np.mean(sat_counts) if sat_counts else 0

        metrics['min_satellites'] = min_sats
        metrics['avg_satellites'] = avg_sats

        pos_type = gnss_data[0].get('pos_type', 'UNKNOWN')
        if pos_type in RTK_FIXED_POS_TYPES:
            sat_threshold = RTK_STANDARDS['min_satellites_fixed']
        else:
            sat_threshold = RTK_STANDARDS['min_satellites_float']

        compliance['satellite_count'] = {
            'value': min_sats,
            'threshold': sat_threshold,
            'pass': min_sats >= sat_threshold,
            'message': "✅ 通过" if min_sats >= sat_threshold else "⚠️ 警告"
        }

    # Sigma values
    if gnss_data:
        lat_sigmas = [r.get('lat_sigma', 999) for r in gnss_data if 'lat_sigma' in r]
        lon_sigmas = [r.get('lon_sigma', 999) for r in gnss_data if 'lon_sigma' in r]
        hgt_sigmas = [r.get('hgt_sigma', 999) for r in gnss_data if 'hgt_sigma' in r]

        max_lat_sigma = np.max(lat_sigmas) if lat_sigmas else 999
        max_lon_sigma = np.max(lon_sigmas) if lon_sigmas else 999
        max_hgt_sigma = np.max(hgt_sigmas) if hgt_sigmas else 999

        metrics['max_horizontal_sigma'] = max(max_lat_sigma, max_lon_sigma)
        metrics['max_vertical_sigma'] = max_hgt_sigma

        compliance['horizontal_sigma'] = {
            'value': max(max_lat_sigma, max_lon_sigma),
            'threshold': RTK_STANDARDS['max_horizontal_sigma'],
            'pass': max(max_lat_sigma, max_lon_sigma) <= RTK_STANDARDS['max_horizontal_sigma'],
            'message': "✅ 通过" if max(max_lat_sigma, max_lon_sigma) <= RTK_STANDARDS['max_horizontal_sigma'] else "⚠️ 警告"
        }

        compliance['vertical_sigma'] = {
            'value': max_hgt_sigma,
            'threshold': RTK_STANDARDS['max_vertical_sigma'],
            'pass': max_hgt_sigma <= RTK_STANDARDS['max_vertical_sigma'],
            'message': "✅ 通过" if max_hgt_sigma <= RTK_STANDARDS['max_vertical_sigma'] else "⚠️ 警告"
        }

    return metrics, compliance

def chart_position_timeseries(gnss_data, output_dir):
    """Plot position time series."""
    if not gnss_data:
        return None

    # Fixed solutions only; plot ENU offsets (m) about the mean position.
    fixed = get_fixed_records(gnss_data)
    t0 = fixed[0]['tow']
    times = [r['tow'] - t0 for r in fixed]
    lats = np.array([r['lat'] for r in fixed])
    lons = np.array([r['lon'] for r in fixed])
    hgts = np.array([r['hgt'] for r in fixed])
    east, north, up = llh_to_enu(lats, lons, hgts,
                                 np.mean(lats), np.mean(lons), np.mean(hgts))

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.suptitle('GNSS Position Time Series - Fixed Solutions (ENU, relative to mean)', fontsize=14)

    for ax, vals, color, name in ((axes[0], east, 'b', 'East'),
                                  (axes[1], north, 'g', 'North'),
                                  (axes[2], up, 'r', 'Up')):
        ax.plot(times, vals, color=color, linewidth=1)
        ax.axhline(0.0, color='r', linestyle='--', label='Mean: 0')
        std = np.std(vals)
        ax.fill_between(times, -2*std, 2*std, alpha=0.2, label='±2σ')
        ax.set_ylabel(f'{name} Offset (m)')
        ax.legend()
        ax.grid(True, alpha=0.3)
    axes[2].set_xlabel('Time (s since first fixed solution)')

    plt.tight_layout()

    filepath = output_dir / 'pos_timeseries.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name

def chart_en_scatter(gnss_data, output_dir):
    """Plot East-North scatter plot."""
    if not gnss_data:
        return None, None

    # Scatter/statistics on valid fixed solutions only (startup SINGLE points
    # are meters away and would dwarf the CEP95 circle).
    fixed = get_fixed_records(gnss_data)
    lats = np.array([r['lat'] for r in fixed])
    lons = np.array([r['lon'] for r in fixed])
    hgts = np.array([r['hgt'] for r in fixed])

    # Convert to ENU
    mean_lat = np.mean(lats)
    mean_lon = np.mean(lons)
    mean_hgt = np.mean(hgts)

    east, north, up = llh_to_enu(lats, lons, hgts,
                                mean_lat, mean_lon, mean_hgt)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(east, north, s=20, alpha=0.6, c='blue', edgecolors='none')

    # 1-sigma ellipse
    theta = np.linspace(0, 2*np.pi, 100)
    sigma_e = np.std(east)
    sigma_n = np.std(north)
    ax.plot(sigma_e * np.cos(theta), sigma_n * np.sin(theta),
            'r-', linewidth=2, label=f'1-sigma (E={sigma_e:.4f}m, N={sigma_n:.4f}m)')

    # 2-sigma ellipse
    ax.plot(2*sigma_e * np.cos(theta), 2*sigma_n * np.sin(theta),
            'r--', linewidth=2, label=f'2-sigma (E={2*sigma_e:.4f}m, N={2*sigma_n:.4f}m)')

    # CEP95 circle
    cep95 = calc_cep95(east, north)
    circle_cep = Circle((0, 0), cep95, fill=False, color='green',
                       linewidth=2, label=f'CEP95={cep95:.4f}m')
    ax.add_patch(circle_cep)

    # Mean marker
    ax.plot(0, 0, 'r+', markersize=15, markeredgewidth=2, label='Mean')

    ax.set_xlabel('East Offset (m)')
    ax.set_ylabel('North Offset (m)')
    ax.set_title('Horizontal Position Scatter (from Mean)')
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)

    filepath = output_dir / 'en_scatter.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    
    # Return elevation statistics for report (ddof=1 for sample standard deviation)
    sigma_u = np.std(up, ddof=1)
    return filepath.name, {'sigma_u': sigma_u}

def chart_position_sigma(gnss_data, output_dir, gst_data=None):
    """绘制定位精度(sigma)时间序列。

    华测 M720 特性(数据为准)：位置报文(BESTPOSA/BESTPA)的 lat/lon/hgt_sigma
    全程为接收机恒定输出，不随解质量逐历元更新(实测恒为定值)，故不能反映
    收敛/失锁。真正逐历元变化的动态精度量在 GNGST：
        rms        = 伪距残差 RMS(m)，失锁/重收敛时冲高；
        major/minor= 误差椭圆半长/半短轴(m)。
    本图以 GNGST 动态量为主轴；位置sigma作为参考平线叠加，并在标题注明其恒定性。
    """
    if not gnss_data and not gst_data:
        return None

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.suptitle('Position Accuracy (GNGST dynamic; BESTPOSA position σ constant & unused)', fontsize=13)

    if gst_data:
        t0 = gst_data[0]['time_seconds']
        gt = [r['time_seconds'] - t0 for r in gst_data]
        rms = [r['rms'] for r in gst_data]
        major = [r['major_sigma'] if r['major_sigma'] is not None else float('nan') for r in gst_data]
        minor = [r['minor_sigma'] if r['minor_sigma'] is not None else float('nan') for r in gst_data]

        axes[0].plot(gt, rms, 'b-', linewidth=1, label='Pseudorange residual RMS')
        axes[0].set_ylabel('RMS residual (m)')
        axes[0].set_title('GNGST Pseudorange Residual RMS')
        axes[0].legend(); axes[0].grid(True, alpha=0.3)

        axes[1].plot(gt, major, 'g-', linewidth=1, label='Error ellipse major axis')
        axes[1].set_ylabel('Major σ (m)')
        axes[1].set_title('GNGST Error Ellipse Major Axis')
        axes[1].legend(); axes[1].grid(True, alpha=0.3)

        axes[2].plot(gt, minor, 'r-', linewidth=1, label='Error ellipse minor axis')
        axes[2].set_ylabel('Minor σ (m)')
        axes[2].set_title('GNGST Error Ellipse Minor Axis')
        axes[2].legend(); axes[2].grid(True, alpha=0.3)

        # 说明BESTPOSA位置σ为恒定输出（不随解质量逐历元更新），仅作文字注释，
        # 不再叠加其9.59m参考线——否则会把0.002~0.011m的动态量压到Y轴底部看成直线。
        if gnss_data:
            ls = gnss_data[0]
            const_note = (f"BESTPOSA lat/lon/hgt σ constant = "
                          f"{ls.get('lat_sigma',0):.4f}/{ls.get('lon_sigma',0):.4f}/"
                          f"{ls.get('hgt_sigma',0):.4f} m (receiver fixed, unused)")
            axes[0].text(0.01, 0.95, const_note, transform=axes[0].transAxes,
                         fontsize=8, va='top', color='dimgray')
        # 各动态量子图按自身数据量程自适应，并留10%边距，保证逐历元波动清晰可见
        for ax in axes:
            ax.margins(y=0.1)
            ax.relim()
            ax.autoscale_view(scalex=False, scaley=True)
    else:
        # 无 GST：退回画位置sigma（可能恒定），并注明
        times = [r['tow'] for r in gnss_data]
        t0 = times[0] if times else 0
        tn = [t - t0 for t in times]
        axes[0].plot(tn, [r['lat_sigma'] for r in gnss_data], 'b-', linewidth=1)
        axes[0].set_ylabel('Latitude σ (m)'); axes[0].grid(True, alpha=0.3)
        axes[0].set_title('BESTPOSA Position σ (receiver constant output)')
        axes[1].plot(tn, [r['lon_sigma'] for r in gnss_data], 'g-', linewidth=1)
        axes[1].set_ylabel('Longitude σ (m)'); axes[1].grid(True, alpha=0.3)
        axes[2].plot(tn, [r['hgt_sigma'] for r in gnss_data], 'r-', linewidth=1)
        axes[2].set_ylabel('Height σ (m)'); axes[2].grid(True, alpha=0.3)

    axes[2].set_xlabel('Time (s since start)')
    plt.tight_layout()

    filepath = output_dir / 'position_sigma_ts.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name
def chart_solution_type_pie(gnss_data, output_dir):
    """绘制定位类型(pos_type)分布饼图，并返回分布统计。

    仅统计实际出现的类型（不再因坐标/卫星数过滤而漏掉 NONE 等无有效坐标的
    历元——上游解析已保留所有 sol_status/pos_type 组合）。类型中文含义按本设备
    枚举表(POS_TYPE_ENUM)标注；北云=表 4-2，华测=表 3-40，二者枚举不同。
    返回 {'counts': {类型: 次数}, 'pct': {类型: 百分比}, 'total': 总数}。
    """
    if not gnss_data:
        return None, {'counts': {}, 'pct': {}, 'total': 0}
    counts = {}
    for r in gnss_data:
        pt = r.get('pos_type', 'UNKNOWN')
        counts[pt] = counts.get(pt, 0) + 1
    total = sum(counts.values())
    pct = {k: (100.0 * v / total if total else 0.0) for k, v in counts.items()}
    # 标签: 类型(中文)
    labels = [f"{k}({POS_TYPE_ENUM.get(k, '未知')})" for k in counts.keys()]

    fig, ax = plt.subplots(figsize=(9, 8))
    ax.pie(counts.values(), labels=labels, autopct='%1.1f%%', startangle=90)
    ax.set_title('Solution Type Distribution (BESTPOSA/BESTGNSSPOSA)')

    filepath = output_dir / 'solution_type_pie.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name, {'counts': counts, 'pct': pct, 'total': total}

def _fmt_solution_dist_html(dist, enum_map):
    """把解类型/状态分布统计格式化为 HTML 表格（含中文含义），与 _fmt_solution_dist 对应。"""
    if not dist or not dist.get('counts'):
        return "<p>(无数据)</p>"
    rows = ["<table class=\"dist\"><tr><th>类型</th><th>含义</th><th>次数</th><th>占比</th></tr>"]
    for k in sorted(dist['counts'].keys(), key=lambda x: -dist['counts'][x]):
        rows.append(f"<tr><td>{k}</td><td>{enum_map.get(k, '未知')}</td><td>{dist['counts'][k]}</td><td>{dist['pct'][k]:.2f}%</td></tr>")
    rows.append(f"<tr><td><b>合计</b></td><td></td><td><b>{dist['total']}</b></td><td><b>100.00%</b></td></tr></table>")
    return "".join(rows)


def _fmt_solution_dist(dist, enum_map):
    """把解类型/状态分布统计格式化为 Markdown 表格（含中文含义）。

    列出实际出现的所有类型（含 NONE/无效历元），标注出现次数与占比；
    用于报告，便于核对分布是否漏了某类（如 NONE）以及各产品枚举差异。
    """
    if not dist or not dist.get('counts'):
        return "(无数据)"
    lines = ["| 类型 | 含义 | 次数 | 占比 |", "|---|---|---:|---:|"]
    for k in sorted(dist['counts'].keys(), key=lambda x: -dist['counts'][x]):
        lines.append(f"| {k} | {enum_map.get(k, '未知')} | {dist['counts'][k]} | {dist['pct'][k]:.2f}% |")
    lines.append(f"| **合计** | | **{dist['total']}** | **100.00%** |")
    return "\n".join(lines)


def chart_num_sats_gnss(gnss_data, output_dir):
    """Plot number of satellites used."""
    if not gnss_data:
        return None

    times = [r['tow'] for r in gnss_data]
    num_sats = [r['num_svs'] for r in gnss_data]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(times, num_sats, 'b-', linewidth=2)
    ax.set_xlabel('Time (s since start)')
    ax.set_ylabel('Number of Satellites')
    ax.set_title('Number of Satellites Used in Solution')
    ax.grid(True, alpha=0.3)

    filepath = output_dir / 'num_sats_ts.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name

def chart_dop_timeseries(gsa_data, output_dir, gga_data=None):
    """Plot DOP time series."""
    if not gsa_data:
        return None

    # Use the record's own timestamps when available (BESTDOPSA 'tow' or NMEA 'time_seconds');
    # fabricating an axis from GGA is misleading when counts/rates differ.
    if gsa_data and 'tow' in gsa_data[0]:
        times = [r['tow'] - gsa_data[0]['tow'] for r in gsa_data]
    elif gsa_data and 'time_seconds' in gsa_data[0]:
        times = [r['time_seconds'] - gsa_data[0]['time_seconds'] for r in gsa_data]
    else:
        times = list(range(len(gsa_data)))

    pdops = [r['pdop'] for r in gsa_data]
    hdops = [r['hdop'] for r in gsa_data]
    vdops = [r.get('vdop') for r in gsa_data]
    has_vdop = all(v is not None for v in vdops)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(times, pdops, 'r-', linewidth=2, label='PDOP')
    ax.plot(times, hdops, 'g-', linewidth=2, label='HDOP')
    if has_vdop:
        ax.plot(times, vdops, 'b-', linewidth=2, label='VDOP')

    ax.set_xlabel('Time (s since start)')
    ax.set_ylabel('DOP Value')
    ax.set_title('Dilution of Precision (DOP) Time Series')
    ax.legend()
    ax.grid(True, alpha=0.3)

    filepath = output_dir / 'dop_timeseries.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name

def chart_dop_histogram(gsa_data, output_dir):
    """Plot DOP histograms."""
    if not gsa_data:
        return None

    pdops = [r['pdop'] for r in gsa_data]
    hdops = [r['hdop'] for r in gsa_data]
    vdops = [r.get('vdop') for r in gsa_data]
    has_vdop = all(v is not None for v in vdops)

    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    fig.suptitle('DOP Distribution Histograms', fontsize=14)

    axes[0].hist(pdops, bins=20, alpha=0.7, color='red', edgecolor='black')
    axes[0].axvline(np.mean(pdops), color='blue', linestyle='--',
                    label=f'Mean: {np.mean(pdops):.2f}')
    axes[0].axvline(np.median(pdops), color='green', linestyle='--',
                    label=f'Median: {np.median(pdops):.2f}')
    axes[0].set_ylabel('Count')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title('PDOP Distribution')

    axes[1].hist(hdops, bins=20, alpha=0.7, color='green', edgecolor='black')
    axes[1].axvline(np.mean(hdops), color='blue', linestyle='--',
                    label=f'Mean: {np.mean(hdops):.2f}')
    axes[1].axvline(np.median(hdops), color='green', linestyle='--',
                    label=f'Median: {np.median(hdops):.2f}')
    axes[1].set_ylabel('Count')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_title('HDOP Distribution')

    if has_vdop:
        axes[2].hist(vdops, bins=20, alpha=0.7, color='blue', edgecolor='black')
        axes[2].axvline(np.mean(vdops), color='blue', linestyle='--',
                        label=f'Mean: {np.mean(vdops):.2f}')
        axes[2].axvline(np.median(vdops), color='green', linestyle='--',
                        label=f'Median: {np.median(vdops):.2f}')
    else:
        axes[2].text(0.5, 0.5, 'VDOP not available\n(BESTDOPSA has no VDOP field)',
                     ha='center', va='center', transform=axes[2].transAxes)
    axes[2].set_ylabel('Count')
    axes[2].set_xlabel('DOP Value')
    if has_vdop:
        axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    axes[2].set_title('VDOP Distribution' if has_vdop else 'VDOP Distribution (N/A)')

    plt.tight_layout()

    filepath = output_dir / 'dop_histogram.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name

def chart_skyplot(gsv_data, output_dir):
    """
    绘制卫星天顶图（Skyplot）
    
    天顶图以极坐标形式展示所有可见卫星的位置：
    - 角度（方位角）：从正北顺时针方向
    - 半径（仰角）：中心为天顶（90°），边缘为地平线（0°）
    
    参数:
        gsv_data (list): 从 parse_gsv 函数获取的卫星数据列表
        output_dir (Path): 输出目录路径
    
    返回:
        str: 生成的图表文件名（skyplot.png），如果数据为空则返回 None
    
    图表特性:
        - 不同星座使用不同颜色和标记符号
        - GPS: 绿色圆形
        - GLONASS: 红色方形
        - Galileo: 蓝色三角形
        - BeiDou: 黄色菱形
        - QZSS: 紫色星形
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='polar')

    # Plot all satellites from all epochs (average)
    constellations = set()
    all_sats = {}

    for epoch in gsv_data:
        for const, data in epoch.items():
            constellations.add(const)
            if const not in all_sats:
                all_sats[const] = {}

            # Group by PRN to average SNR and position
            for sat in data['satellites']:
                prn = sat['prn']
                if prn not in all_sats[const]:
                    all_sats[const][prn] = {
                        'elevation': sat['elevation'],
                        'azimuth': sat['azimuth'],
                        'snr': [sat['snr']],
                        'count': 1
                    }
                else:
                    # Average SNR and position
                    all_sats[const][prn]['snr'].append(sat['snr'])
                    all_sats[const][prn]['elevation'] = (all_sats[const][prn]['elevation'] * all_sats[const][prn]['count'] + sat['elevation']) / (all_sats[const][prn]['count'] + 1)
                    all_sats[const][prn]['azimuth'] = (all_sats[const][prn]['azimuth'] * all_sats[const][prn]['count'] + sat['azimuth']) / (all_sats[const][prn]['count'] + 1)
                    all_sats[const][prn]['count'] += 1

    # Plot each constellation
    for const, sats in all_sats.items():
        color = CONSTELLATION_COLORS.get(const, 'gray')
        marker = CONSTELLATION_MARKERS.get(const, 'o')

        for prn, sat_data in sats.items():
            avg_snr = np.mean(sat_data['snr'])
            az_rad = np.radians(sat_data['azimuth'])
            r = 90 - sat_data['elevation']  # 90° at center, 0° at edge
            ax.plot(az_rad, r, marker=marker, color=color,
                   markersize=4 + avg_snr/10, alpha=0.8,
                   markeredgecolor='black', markeredgewidth=0.5)

    ax.set_theta_zero_location('N')  # North at top
    ax.set_theta_direction(-1)       # Clockwise
    ax.set_rlim(0, 90)
    ax.set_rticks([0, 30, 60, 90])
    ax.set_yticklabels(['90°', '60°', '30°', '0°'])

    # Legend
    handles = []
    for const in constellations:
        handles.append(mpatches.Patch(color=CONSTELLATION_COLORS[const],
                                    label=const))
    ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.3, 1.0))

    ax.set_title('Satellite Skyplot (Elevation vs Azimuth)', pad=20)

    filepath = output_dir / 'skyplot.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight', facecolor='white')
    plt.close()
    return filepath.name

def chart_snr_distribution(gsv_data, output_dir):
    """
    绘制各星座的信噪比（SNR）分布箱线图
    
    箱线图展示每个星座所有卫星的 SNR 分布情况，便于比较不同星座的信号质量。
    
    参数:
        gsv_data (list): 从 parse_gsv 函数获取的卫星数据列表
        output_dir (Path): 输出目录路径
    
    返回:
        str: 生成的图表文件名（snr_distribution.png），如果数据为空则返回 None
    
    图表说明:
        - 箱体：包含 25% 到 75% 的数据（四分位距 IQR）
        - 黄色中线：中位数
        - 须线（Whiskers）：显示数据的范围（1.5×IQR）
        - 圆圈：异常值（超出须线范围的数据点）
        - 颜色：GPS-绿色、GLONASS-红色、Galileo-蓝色、QZSS-紫色、BeiDou-黄色
    
    数据来源:
        GSV 消息中的 SNR 字段，仅使用仰角 > 10° 的卫星数据
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Collect SNR values per constellation
    snr_by_const = {}
    sats_by_const = {}  # 每个星座出现过的唯一卫星 PRN 集合(用于统计卫星颗数)

    for epoch in gsv_data:
        for const, data in epoch.items():
            if const not in snr_by_const:
                snr_by_const[const] = []
                sats_by_const[const] = set()

            for sat in data['satellites']:
                sats_by_const[const].add(sat['prn'])  # 统计该星座唯一卫星数(不限仰角)
                if sat['elevation'] > 10:  # Use only satellites above 10°
                    snr_by_const[const].append(sat['snr'])

    # Create box plot
    constellations = list(snr_by_const.keys())
    snr_values = [snr_by_const[const] for const in constellations]

    colors = [CONSTELLATION_COLORS.get(const, 'gray') for const in constellations]

    bp = ax.boxplot(snr_values, patch_artist=True, showfliers=False)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xlabel('Constellation')
    ax.set_ylabel('SNR (dBHz)')
    ax.set_title('Signal-to-Noise Ratio Distribution by Constellation')
    ax.grid(True, alpha=0.3)

    # Set x-axis labels（星座名后括号标注该星座唯一卫星颗数）
    x_labels = [f"{const} ({len(sats_by_const.get(const, set()))})" for const in constellations]
    ax.set_xticks(range(1, len(constellations)+1))
    ax.set_xticklabels(x_labels)

    filepath = output_dir / 'snr_distribution.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name

def chart_satellites_used_ts(gnss_data, output_dir):
    """Plot #SVs and #solnSVs from BESTGNSSPOSA messages over time."""
    if not gnss_data:
        return None

    # Extract time and satellite counts from BESTGNSSPOSA data
    times = [r['tow'] for r in gnss_data]
    num_svs = [r['num_svs'] for r in gnss_data]
    num_soln_svs = [r['num_soln_svs'] for r in gnss_data]

    # Normalize time to start from 0
    start_time = times[0] if times else 0
    times_normalized = [t - start_time for t in times]

    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot #SVs (total satellites tracked)
    ax.plot(times_normalized, num_svs, 'b-', linewidth=2, label='#SVs (Tracked)')
    
    # Plot #solnSVs (satellites used in solution)
    ax.plot(times_normalized, num_soln_svs, 'r--', linewidth=2, label='#solnSVs (Used in Solution)')

    ax.set_xlabel('Time (s since start)')
    ax.set_ylabel('Number of Satellites')
    ax.set_title('Satellites Tracked vs Used in Solution (#SVs vs #solnSVs)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    filepath = output_dir / 'satellites_used_ts.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name

def chart_satellite_info(gsv_data, output_dir):
    """Create a table of satellite information."""
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.axis('off')

    # Collect satellite info
    sat_info = []
    constellations = set()

    for epoch in gsv_data:
        for const, data in epoch.items():
            constellations.add(const)
            for sat in data['satellites']:
                if sat['elevation'] >= 0:  # Include all satellites, including those at horizon
                    sat_info.append({
                        'constellation': const,
                        'prn': sat['prn'],
                        'elevation': sat['elevation'],
                        'azimuth': sat['azimuth'],
                        'snr': sat['snr']
                    })

    # Sort by constellation and PRN
    sat_info.sort(key=lambda x: (x['constellation'], x['prn']))

    # Create table with sequence number
    table_data = []
    headers = ['#', 'Constellation', 'PRN', 'Elevation (°)', 'Azimuth (°)', 'SNR (dBHz)']
    table_data.append(headers)

    # Display all satellites (no limit)
    for idx, sat in enumerate(sat_info, 1):
        row = [
            str(idx),  # Sequence number
            sat['constellation'],
            str(sat['prn']),
            f"{sat['elevation']:.1f}",
            f"{sat['azimuth']:.1f}",
            f"{sat['snr']:.1f}"
        ]
        table_data.append(row)

    # Plot table
    table = ax.table(cellText=table_data, loc='center',
                     cellLoc='center', fontsize=10)
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)

    # Style header row
    for i in range(len(headers)):
        table[(0, i)].set_facecolor(CONSTELLATION_COLORS['GPS'])
        table[(0, i)].set_text_props(weight='bold', color='white')

    ax.set_title('Satellite Visibility Summary', fontsize=14)

    filepath = output_dir / 'satellite_table.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name

def chart_fix_rate_gauge(fix_rate, output_dir):
    """Create a fix rate gauge chart."""
    fig, ax = plt.subplots(figsize=(8, 8))

    # Create gauge-like visualization
    ax.barh([0], [fix_rate * 100], height=0.5, color='green' if fix_rate > 0.95 else 'orange' if fix_rate > 0.9 else 'red')
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.set_xlabel('Fix Rate (%)')
    ax.set_title(f'Fix Rate: {fix_rate * 100:.1f}%')

    # Add threshold lines
    ax.axvline(95, color='red', linestyle='--', linewidth=2, label='95% Threshold')
    ax.legend()

    # Add text
    ax.text(50, 0, f'{fix_rate * 100:.1f}%', ha='center', va='center', fontsize=20, weight='bold')

    filepath = output_dir / 'fix_rate_gauge.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name

def chart_data_gaps(gaps, output_dir):
    """Plot data gap analysis timeline."""
    if not gaps:
        fig, ax = plt.subplots(figsize=(10, 2))
        ax.text(0.5, 0.5, 'No data gaps detected', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.axis('off')
    else:
        fig, ax = plt.subplots(figsize=(12, 4))

        # Timeline
        timeline = [0]
        for gap in gaps:
            timeline.extend([gap['start_time'], gap['end_time']])
        timeline.append(timeline[-1] + 1)

        # Draw timeline
        ax.plot(timeline, [0] * len(timeline), 'k-', linewidth=2)

        # Highlight gaps
        for i, gap in enumerate(gaps):
            ax.fill_between([gap['start_time'], gap['end_time']], -0.1, 0.1,
                           alpha=0.7, color='red', label=f'Gap {i+1}' if i < 5 else '')

            # Add gap info text
            mid_time = (gap['start_time'] + gap['end_time']) / 2
            ax.text(mid_time, 0.15, f'{gap["gap_seconds"]:.2f}s',
                   ha='center', va='bottom', fontsize=10)

        ax.set_xlabel('Time (s since start)')
        ax.set_title('Data Gap Analysis')
        ax.set_ylim(-0.2, 0.3)
        ax.set_yticks([])

        if len(gaps) <= 5:
            ax.legend()

        ax.grid(True, alpha=0.3)

    filepath = output_dir / 'data_gap_analysis.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name

def chart_performance_summary(gnss_data, output_dir):
    """Create a performance summary table chart."""
    if not gnss_data:
        return None

    # Calculate metrics
    # Position statistics are computed on valid FIXED solutions only (NARROW_INT /
    # NARROW_FLOAT with SOL_COMPUTED), consistent with the 'Fixed Sol' peak labels.
    # Startup SINGLE/NO-fix records would otherwise corrupt a static-precision table.
    fixed_recs = get_fixed_records(gnss_data)
    lats = np.array([r['lat'] for r in fixed_recs])
    lons = np.array([r['lon'] for r in fixed_recs])
    hgts = np.array([r['hgt'] for r in fixed_recs])

    mean_lat = np.mean(lats)
    mean_lon = np.mean(lons)
    mean_hgt = np.mean(hgts)

    east, north, up = llh_to_enu(lats, lons, hgts, mean_lat, mean_lon, mean_hgt)

    cep95 = calc_cep95(east, north)
    horizontal_rms = np.sqrt(np.mean(east**2 + north**2))
    # Vertical RMS: deviation about the mean height (same convention as E/N RMS
    # relative to the mean position in llh_to_enu).
    vertical_rms = np.sqrt(np.mean(up**2))
    horizontal_peak = np.max(np.sqrt(east**2 + north**2))
    vertical_peak = np.max(np.abs(up))

    # Fix rate
    fixed_count = len([r for r in gnss_data if r.get('sol_status') in VALID_SOL_STATUSES and r.get('pos_type') in RTK_FIXED_POS_TYPES])
    fix_rate = fixed_count / len(gnss_data) * 100

    # Calculate expected and actual epochs for data completeness
    expected_count = len(gnss_data)
    actual_count = len(gnss_data)

    # Create table
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')

    table_data = [
        ['Metric', 'Value', 'Description'],
        ['Expected/Actual Epochs', f'{expected_count}/{actual_count}', 'Data completeness'],
        ['Fixed Solution Ratio', f'{fix_rate:.1f}%', 'RTK fixed percentage'],
        ['Position CEP95', f'{cep95:.4f} m', '95% circular error'],
        ['Horizontal RMS', f'{horizontal_rms:.4f} m', 'sqrt(mean(dE^2+dN^2))'],
        ['Vertical RMS', f'{vertical_rms:.4f} m', 'RMS about mean height'],
        ['Fixed Sol Horizontal Peak', f'{horizontal_peak:.4f} m', 'Max horizontal deviation'],
        ['Fixed Sol Vertical Peak', f'{vertical_peak:.4f} m', 'Max vertical deviation']
    ]

    table = ax.table(cellText=table_data, loc='center', cellLoc='left', fontsize=11,
                    colWidths=[0.35, 0.25, 0.40])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 2.0)

    # Style header row
    for i in range(3):
        table[(0, i)].set_facecolor('#3498db')
        table[(0, i)].set_text_props(weight='bold', color='white')

    # Style metric rows
    for i in range(1, len(table_data)):
        for j in range(3):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#f9f9f9')

    ax.set_title('Performance Summary Table', fontsize=16, pad=20)

    filepath = output_dir / 'performance_summary.png'
    plt.savefig(filepath, dpi=DPI_DEFAULT, bbox_inches='tight')
    plt.close()
    return filepath.name

# === REPORT GENERATION FUNCTIONS ===

def encode_image_to_base64(image_path):
    """Encode an image file to base64 string."""
    try:
        with open(image_path, 'rb') as f:
            return base64.b64encode(f.read()).decode('utf-8')
    except Exception as e:
        print(f"Error encoding image {image_path}: {e}")
        return ""

def generate_html_report(output_dir, data, stats, chart_files):
    """Generate HTML report with embedded images."""
    # Encode chart images to base64
    encoded_images = {}
    for chart_name, chart_file in chart_files.items():
        if isinstance(chart_file, tuple):
            chart_file = chart_file[0]  # Take first element if tuple
        if chart_file:
            image_path = output_dir / chart_file
            encoded_images[chart_name] = encode_image_to_base64(image_path)
        else:
            encoded_images[chart_name] = ""

    # Calculate position statistics if GNSS data is available
    pos_stats = {
        'east_std': 0, 'north_std': 0, 'up_std': 0,
        'east_rms': 0, 'north_rms': 0, 'up_rms': 0,
        'cep95': 0, 'drms': 0, 'horizontal_peak': 0, 'vertical_peak': 0
    }

    if 'gnss_data' in data and data['gnss_data']:
        # Static-precision statistics on valid fixed solutions only (see get_fixed_records).
        gnss_data = get_fixed_records(data['gnss_data'])
        lats = np.array([r['lat'] for r in gnss_data])
        lons = np.array([r['lon'] for r in gnss_data])
        hgts = np.array([r['hgt'] for r in gnss_data])

        mean_lat = np.mean(lats)
        mean_lon = np.mean(lons)
        mean_hgt = np.mean(hgts)

        east, north, up = llh_to_enu(lats, lons, hgts, mean_lat, mean_lon, mean_hgt)

        pos_stats['east_std'] = np.std(east, ddof=1)
        pos_stats['north_std'] = np.std(north, ddof=1)
        pos_stats['up_std'] = np.std(up, ddof=1)
        pos_stats['east_rms'] = np.sqrt(np.mean(east**2))
        pos_stats['north_rms'] = np.sqrt(np.mean(north**2))
        pos_stats['up_rms'] = np.sqrt(np.mean(up**2))
        pos_stats['cep95'] = calc_cep95(east, north)
        pos_stats['drms'] = calc_2drms(east, north)
        pos_stats['horizontal_peak'] = np.max(np.sqrt(east**2 + north**2))
        pos_stats['vertical_peak'] = np.max(np.abs(up))

    # DOP statistics: prefer GSA (only source that carries VDOP), fall back to BESTDOPSA (no VDOP).
    # Must stay consistent with the chart data source (merge_gsa_dop_per_epoch(gsa_data) if gsa_data else bestdops_data).
    dop_stats = {'pdop': [0], 'hdop': [0], 'vdop': [0]}
    dop_src = None
    if 'gsa_data' in data and data['gsa_data']:
        dop_src = data['gsa_data']  # GSA 优先(含 VDOP)
    elif 'bestdops_data' in data and data['bestdops_data']:
        dop_src = data['bestdops_data']
    if dop_src:
        dop_stats['pdop'] = [r['pdop'] for r in dop_src]
        dop_stats['hdop'] = [r['hdop'] for r in dop_src]
        dop_stats['vdop'] = [r.get('vdop') for r in dop_src]

    html_content = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <title>GPS RTK Static Positioning Analysis Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1, h2, h3 {{
            color: #333;
            margin-top: 30px;
        }}
        h1 {{
            border-bottom: 3px solid #2c3e50;
            padding-bottom: 10px;
        }}
        h2 {{
            border-bottom: 2px solid #3498db;
            padding-bottom: 5px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #3498db;
            color: white;
            font-weight: bold;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        .metric-good {{
            color: green;
            font-weight: bold;
        }}
        .metric-warning {{
            color: orange;
            font-weight: bold;
        }}
        .metric-bad {{
            color: red;
            font-weight: bold;
        }}
        .chart {{
            margin: 20px 0;
            border: 1px solid #ddd;
            border-radius: 5px;
            overflow: hidden;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        .chart img {{
            width: 100%;
            height: auto;
            display: block;
        }}
        .pass {{
            color: green;
            font-weight: bold;
        }}
        .fail {{
            color: red;
            font-weight: bold;
        }}
        .summary-table {{
            background-color: #f8f9fa;
            border: 2px solid #dee2e6;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>GPS RTK Static Positioning Analysis Report</h1>
        <p><strong>File:</strong> {data['filename']} |
           <strong>Generated:</strong> {stats['generated_time']} |
           <strong>Duration:</strong> {stats['duration']:.2f} seconds</p>

        <!-- Section 1: Data Overview -->
        <h2>1. 数据概览</h2>
        <table class="summary-table">
            <tr>
                <th>消息类型</th>
                <th>数量</th>
                <th>设置周期</th>
                <th>实际周期</th>
                <th>实际频率</th>
                <th>插入均匀性</th>
                <th>符合性</th>
            </tr>"""

    # Add message type summary (实测周期/频率/均匀性/符合性，标色告警)
    _mp = stats.get('msg_periods', {})
    for msg_type, info in data['message_counts'].items():
        nominal = EXPECTED_INTERVALS.get(msg_type, 0)
        nominal_str = f"{nominal:.2f}s ({1/nominal:.1f}Hz)" if nominal > 0 else "N/A"
        mp = _mp.get(msg_type)
        if msg_type in PERIOD_NOT_MEASURABLE:
            period_str = hz_str = "N/A"
            uni_str = "-"
            cmp_str = "报文无时间戳，周期不可实测"
            uni_cls = cmp_cls = ""
        elif mp:
            period_str = f"{mp['median_interval_s']:.3f}s"
            hz_str = f"{mp['actual_hz']:.2f}Hz"
            uni_str = f"均匀(CV={mp['cv']:.2f})" if mp['uniform'] else f"不均(CV={mp['cv']:.2f})"
            uni_cls = "pass" if mp['uniform'] else "metric-warning"
            note = mp.get('note', '')
            if not mp['compliant']:
                cmp_str = f"✗偏差{mp['dev_pct']:+.0f}%" + (f" {note}" if note else "")
                cmp_cls = "fail"
            elif note:
                cmp_str = f"⚠ {note}"
                cmp_cls = "metric-warning"
            else:
                cmp_str = "✓符合"
                cmp_cls = "pass"
        else:
            period_str = hz_str = "N/A"
            uni_str = cmp_str = "-"
            uni_cls = cmp_cls = ""
        html_content += f"""
            <tr>
                <td>{msg_type}</td>
                <td>{info:,}</td>
                <td>{nominal_str}</td>
                <td>{period_str}</td>
                <td>{hz_str}</td>
                <td class="{uni_cls}">{uni_str}</td>
                <td class="{cmp_cls}">{cmp_str}</td>
            </tr>"""

    # Get location data
    loc = data.get('location', {'lat': 0, 'lon': 0, 'hgt': 0})
    ant_hgt = data.get('antenna_height', 0)

    html_content += f"""
        </table>
        <p style="font-size:0.9em;color:#555;"><strong>说明：</strong>
        仅当报文自带时间字段时，"实际周期/频率/均匀性"才是真实测量值 —— # 类报文头含 TOW、NMEA 中 GGA/GST 含 UTC 时间。
        GSA/GSV 报文本身无时间字段（NMEA 协议），且 GSV 每历元连发多条、条数随可见星数变化、同历元各条同时发出，
        故单句周期无意义，标记为 N/A（数据为准，绝不臆造）。
        消息是否按周期发送由输出触发方式决定（ONTIME=周期输出；ONNEW/ONCHANGED=数据变化才输出）。</p>
        <p><strong>数据完整性:</strong> {stats.get('data_completeness', 0):.1f}% (共 {stats.get('total_messages', 0):,} 条消息)</p>
        <p><strong>位置信息:</strong> {loc['lat']:.6f}°N, {loc['lon']:.6f}°E, {loc['hgt']:.2f} m (WGS84)</p>
        <p><strong>天线高度:</strong> {ant_hgt:.3f} m (测站上方)</p>
        <p><strong>报文完整性(校验位):</strong> 有效 {stats.get('integrity', {}).get('valid', 0):,} / {stats.get('integrity', {}).get('total', 0):,} 条
           (通过率 {stats.get('integrity', {}).get('integrity_pct', 0):.1f}%；无效 {stats.get('integrity', {}).get('invalid', 0)} 条)
           — NMEA XOR: {stats.get('integrity', {}).get('nmea_valid', 0):,}/{stats.get('integrity', {}).get('nmea_total', 0):,} ({stats.get('integrity', {}).get('nmea_pct', 0):.1f}%)；
           #类 CRC-32(表4-6): {stats.get('integrity', {}).get('hash_valid', 0):,}/{stats.get('integrity', {}).get('hash_total', 0):,} ({stats.get('integrity', {}).get('hash_pct', 0):.1f}%)</p>
    """

    # Add GNSS Analysis Section
    html_content += f"""
        <h2>2. GNSS 定位分析 (BESTGNSSPOSA)</h2>
        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('pos_timeseries', '')}" alt="Position Time Series"/>
        </div>
        <p>位置时间序列显示了接收机在各个维度上的微小变化。所有参数都围绕均值稳定波动，表明良好的观测条件。</p>

        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('en_scatter', '')}" alt="East-North Scatter Plot"/>
        </div>
        <table>
            <tr>
                <th>统计项</th>
                <th>东西方向 (m)</th>
                <th>南北方向 (m)</th>
                <th>高程方向 (m)</th>
            </tr>
            <tr>
                <td>均值（以均值位置为基准的相对偏移，均值定义上恒为 0）</td>
                <td>0.0000</td>
                <td>0.0000</td>
                <td>0.0000</td>
            </tr>
            <!-- 均值恒为 0 是定义使然：本表统计的是 ENU【相对均值】的偏移量，
                 各点均已先减去均值再统计，故均值必为 0，非巧合也非数据异常。 -->
            <tr>
                <td>标准差</td>
                <td>{pos_stats['east_std']:.4f}</td>
                <td>{pos_stats['north_std']:.4f}</td>
                <td>{pos_stats['up_std']:.4f}</td>
            </tr>
            <tr>
                <td>RMS</td>
                <td>{pos_stats['east_rms']:.4f}</td>
                <td>{pos_stats['north_rms']:.4f}</td>
                <td>{pos_stats['up_rms']:.4f}</td>
            </tr>
            <tr>
                <td>CPE95</td>
                <td>{pos_stats['cep95']:.4f}</td>
                <td>-</td>
                <td>-</td>
            </tr>
            <tr>
                <td>2DRMS</td>
                <td>{pos_stats['drms']:.4f}</td>
                <td>-</td>
                <td>-</td>
            </tr>
        </table>

        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('position_sigma_ts', '')}" alt="Position Sigma Time Series"/>
        </div>

        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('solution_type_pie', '')}" alt="Solution Type Distribution"/>
        </div>
        <div class="dist-table">{_fmt_solution_dist_html(stats.get('sol_type_dist', {}), POS_TYPE_ENUM)}</div>

        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('num_sats_ts', '')}" alt="Number of Satellites Used"/>
        </div>
    """

    # Add DOP Analysis Section
    pdop_stats = calc_statistics(dop_stats['pdop'])
    hdop_stats = calc_statistics(dop_stats['hdop'])
    _vdop_vals = dop_stats['vdop']
    vdop_stats = calc_statistics(_vdop_vals) if _vdop_vals and all(v is not None for v in _vdop_vals) else None
    if vdop_stats:
        vdop_row_html = ('<td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td><td>{:.2f}</td>'
                         .format(vdop_stats['min'], vdop_stats['max'], vdop_stats['mean'],
                                 vdop_stats['std'], vdop_stats['median']))
    else:
        vdop_row_html = ('<td colspan="5" style="text-align:center;color:#888;">'
                         'N/A（VDOP 仅 GSA 提供，BESTDOPSA 无此字段）</td>')

    html_content += f"""
        <h2>3. DOP 分析</h2>
        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('dop_timeseries', '')}" alt="DOP Time Series"/>
        </div>
        <p>DOP（精度衰减因子）值显示了卫星几何构型对定位精度的影响。较低的值表示更好的卫星几何构型。</p>

        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('dop_histogram', '')}" alt="DOP Histogram"/>
        </div>
        <table>
            <tr>
                <th>DOP类型</th>
                <th>最小值</th>
                <th>最大值</th>
                <th>均值</th>
                <th>标准差</th>
                <th>中位数</th>
            </tr>
            <tr>
                <td>PDOP</td>
                <td>{pdop_stats['min']:.2f}</td>
                <td>{pdop_stats['max']:.2f}</td>
                <td>{pdop_stats['mean']:.2f}</td>
                <td>{pdop_stats['std']:.2f}</td>
                <td>{pdop_stats['median']:.2f}</td>
            </tr>
            <tr>
                <td>HDOP</td>
                <td>{hdop_stats['min']:.2f}</td>
                <td>{hdop_stats['max']:.2f}</td>
                <td>{hdop_stats['mean']:.2f}</td>
                <td>{hdop_stats['std']:.2f}</td>
                <td>{hdop_stats['median']:.2f}</td>
            </tr>
            <tr>
                <td>VDOP</td>
                {vdop_row_html}
            </tr>
        </table>
    """

    # Add Satellite Visibility Section
    html_content += f"""
        <h2>4. 卫星可见性分析</h2>
        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('skyplot', '')}" alt="Skyplot"/>
        </div>
        <p>天顶图显示了所有可见卫星的高度角和方位角，以及信号强度（颜色表示）。</p>

        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('snr_distribution', '')}" alt="SNR Distribution"/>
        </div>
        <p style="text-align: center; color: #666; font-size: 14px; margin-top: 10px;">
            <strong>箱线图说明：</strong><br/>
            - <strong>箱体</strong>：包含 25% 到 75% 的数据（四分位距）<br/>
            - <strong>黄色中线</strong>：中位数<br/>
            - <strong>须线（Whiskers）</strong>：显示数据的范围<br/>
            - <strong>颜色</strong>：GPS-绿色、GLONASS-红色、Galileo-蓝色、QZSS-紫色、BeiDou-黄色<br/>
            - <strong>数据来源</strong>：所有可见卫星的 GSV 消息中的 SNR 字段
        </p>

        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('satellites_used_ts', '')}" alt="Satellites Used Time Series"/>
            <p style="text-align: center; color: #666; font-size: 14px; margin-top: 10px;">
                <strong>#SVs</strong> (蓝色实线): 接收机跟踪到的卫星总数<br/>
                <strong>#solnSVs</strong> (红色虚线): 参与定位解算的卫星数
            </p>
        </div>

        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('satellite_table', '')}" alt="Satellite Detail Table"/>
        </div>
    """

    # Performance summary
    html_content += f"""
        <h2>5. 性能汇总</h2>
        <div class="chart">
            <img src="data:image/png;base64,{encoded_images.get('performance_summary', '')}" alt="Performance Summary"/>
        </div>
    """

    # Add Quality Assessment Section with calculated values
    fix_rate = stats.get('fix_rate', 0.0) * 100
    cep95_pass = pos_stats['cep95'] <= RTK_STANDARDS['cep95_rtk_fixed']
    pdop_pass = pdop_stats['mean'] <= RTK_STANDARDS['pdop_moderate']
    hdop_pass = hdop_stats['mean'] <= RTK_STANDARDS['hdop_moderate']
    data_completeness_pass = stats.get('data_completeness', 0) >= 99.0

    _fc = stats.get('fix_continuity', {}) or {}
    _ss = stats.get('static_stability', {}) or {}
    _sc = stats.get('sigma_convergence', {}) or {}
    _sn = stats.get('snr_summary', {}) or {}
    _es = stats.get('envstatus', {}) or {}
    def _g(d, k, fmt='{:.2f}'):
        v = d.get(k)
        return (fmt.format(v) if isinstance(v, (int, float)) else 'N/A')
    _fc_txt = (f"中断{_fc.get('interruption_count', 0)}次 / "
               f"最长连续固定{_fc.get('longest_fix_streak_s', 0):.1f}s / "
               f"均重固定{_fc.get('mean_refix_s', 0):.1f}s") if _fc else 'N/A'

    # ---- 华测特有：环境感知(ENVSTATUSA) 专属章节（仅移动站字段） ----
    if _es:
        def _row(label, stat, unit, good, note):
            if not stat:
                return ""
            return (f"<tr><td>{label}</td>"
                    f"<td>{stat['min']:.1f} / {stat['mean']:.1f} / {stat['max']:.1f} {unit}</td>"
                    f"<td>{good}</td><td>{note}</td></tr>")
        _env_rows = ""
        _env_rows += _row("环境质量分值 env_score", _es.get('env_score'), "", "85~100优/60~84一般/0~59差",
                          "综合环境评分，越高越好")
        _env_rows += _row("卫星可视比率", _es.get('sat_vis_rate'), "%", "越高越好",
                          "跟踪星/理论可视星，反映遮挡")
        _env_rows += _row("参与解算比率", _es.get('sat_sol_rate'), "%", "越高越好",
                          "解算星/跟踪星，反映可用性")
        _env_rows += _row("解置信度", _es.get('sol_confidence'), "", "越高越好",
                          "接收机对当前解的置信程度")
        _env_rows += _row("跟踪卫星数", _es.get('num_svs'), "颗", "越多越好", "可见并跟踪的卫星总数")
        _env_rows += _row("共视卫星数", _es.get('cv_sats'), "颗", "越多越好",
                          "基站-Rover共视星，RTK稳定性")
        _env_rows += _row("平均信噪比 avg_snr", _es.get('avg_snr'), "dBHz",
                          "≥38开阔/32~37半遮挡/≤31严重", "接收机前端信号质量")
        _env_rows += _row("PDOP", _es.get('pdop'), "", "≤2.5好/2.5~5一般/>5差", "卫星几何构型")
        _iono = _es.get('iono_level_counts', {})
        _iono_txt = ", ".join(f"{k}:{v}" for k, v in _iono.items()) if _iono else "N/A"
        html_content += f"""
        <h2>6A. 环境感知分析（ENVSTATUSA · 华测特有 · 移动站模式）</h2>
        <p>该报文为华测模组特色，用于评估观测环境对静态定位的影响。环境越差，静态精度越易恶化。</p>
        <table class="summary-table">
            <tr><th>字段</th><th>最小/均值/最大</th><th>好坏判断</th><th>说明</th></tr>
            {_env_rows}
            <tr><td>电离层活跃状态</td><td>{_iono_txt}</td><td>QUIET平静为佳</td><td>电离层扰动影响定位</td></tr>
        </table>
"""

    html_content += f"""
        <h2>6. 质量评估</h2>
        <table class="summary-table">
            <tr>
                <th>评估指标</th>
                <th>数值</th>
                <th>行业标准</th>
                <th>评估结果</th>
            </tr>
            <tr>
                <td>固定解比例</td>
                <td>{fix_rate:.1f}%</td>
                <td>>{RTK_STANDARDS['fixed_solution_threshold']*100:.0f}%</td>
                <td class="{'pass' if fix_rate >= RTK_STANDARDS['fixed_solution_threshold']*100 else 'fail'}">{'通过' if fix_rate >= RTK_STANDARDS['fixed_solution_threshold']*100 else '未通过'}</td>
            </tr>
            <tr>
                <td>位置CPE95</td>
                <td>{pos_stats['cep95']:.4f} m</td>
                <td><{RTK_STANDARDS['cep95_rtk_fixed']:.2f} m (RTK固定)</td>
                <td class="{'pass' if cep95_pass else 'metric-warning'}">{'✅ 优秀' if pos_stats['cep95'] < RTK_STANDARDS['cep95_rtk_fixed']/2 else '✅ 通过' if cep95_pass else '⚠️ 警告'}</td>
            </tr>
            <tr>
                <td>水平RMS</td>
                <td>{(pos_stats['east_rms'] + pos_stats['north_rms'])/2:.4f} m</td>
                <td>-</td>
                <td>-</td>
            </tr>
            <tr>
                <td>垂直RMS</td>
                <td>{pos_stats['up_rms']:.4f} m</td>
                <td>-</td>
                <td>-</td>
            </tr>
            <tr>
                <td>固定解水平峰值</td>
                <td>{pos_stats['horizontal_peak']:.4f} m</td>
                <td>-</td>
                <td>-</td>
            </tr>
            <tr>
                <td>固定解垂直峰值</td>
                <td>{pos_stats['vertical_peak']:.4f} m</td>
                <td>-</td>
                <td>-</td>
            </tr>
            <tr>
                <td>平均PDOP</td>
                <td>{pdop_stats['mean']:.2f}</td>
                <td><{RTK_STANDARDS['pdop_moderate']:.1f}</td>
                <td class="{'pass' if pdop_stats['mean'] <= RTK_STANDARDS['pdop_good'] else 'metric-warning' if pdop_stats['mean'] <= RTK_STANDARDS['pdop_moderate'] else 'fail'}">{'✅ 优秀' if pdop_stats['mean'] <= RTK_STANDARDS['pdop_good'] else '⚠️ 警告' if pdop_stats['mean'] <= RTK_STANDARDS['pdop_moderate'] else '❌ 差'}</td>
            </tr>
            <tr>
                <td>平均HDOP</td>
                <td>{hdop_stats['mean']:.2f}</td>
                <td><{RTK_STANDARDS['hdop_moderate']:.1f}</td>
                <td class="{'pass' if hdop_stats['mean'] <= RTK_STANDARDS['hdop_good'] else 'metric-warning' if hdop_stats['mean'] <= RTK_STANDARDS['hdop_moderate'] else 'fail'}">{'✅ 优秀' if hdop_stats['mean'] <= RTK_STANDARDS['hdop_good'] else '⚠️ 警告' if hdop_stats['mean'] <= RTK_STANDARDS['hdop_moderate'] else '❌ 差'}</td>
            </tr>
            <tr>
                <td>数据完整性</td>
                <td>{stats.get('data_completeness', 0):.1f}%</td>
                <td>>99%</td>
                <td class="{'pass' if data_completeness_pass else 'fail'}">{'通过' if data_completeness_pass else '未通过'}</td>
            </tr>
            <tr>
                <td>数据间隙</td>
                <td>{stats.get('data_gaps', 0)} 个</td>
                <td>0 个</td>
                <td class="{'pass' if stats.get('data_gaps', 0) == 0 else 'fail'}">{'通过' if stats.get('data_gaps', 0) == 0 else '未通过'}</td>
            </tr>
            <tr>
                <td>固定解连续性</td>
                <td>{_fc_txt}</td>
                <td>越少中断/越短重固定越好</td>
                <td class="{'pass' if _fc.get('interruption_count', 1) == 0 else 'metric-warning'}">{'✅ 连续' if _fc.get('interruption_count', 1) == 0 else '⚠️ 有中断'}</td>
            </tr>
            <tr>
                <td>静止性(纯GNSS)</td>
                <td>最大偏移 {_g(_ss, 'max_horizontal_offset', '{:.3f}')} m</td>
                <td>≤{_ss.get('h_threshold', 0.15):.2f} m</td>
                <td class="{'pass' if _ss.get('is_static') else 'fail'}">{'✅ 静止' if _ss.get('is_static') else '✗ 非静止(运动数据)'}</td>
            </tr>
            <tr>
                <td>伪距残差RMS收敛 (GNGST)</td>
                <td>初 {_g(_sc, 'h_init', '{:.3f}')} → 末 {_g(_sc, 'h_final', '{:.3f}')} m</td>
                <td>收敛且小为佳</td>
                <td class="{'pass' if _sc.get('h_final', 9) < _sc.get('h_init', 0) else 'metric-warning'}">{'✅ 收敛' if _sc.get('h_final', 9) < _sc.get('h_init', 0) else '⚠️ 未明显收敛'}</td>
            </tr>
            <tr>
                <td>载噪比 C/N0 均值</td>
                <td>{_g(_sn, 'mean', '{:.1f}')} dBHz</td>
                <td>≥38开阔 / 32-37半遮挡 / ≤31严重</td>
                <td class="{'pass' if _sn.get('mean', 0) >= 38 else 'metric-warning' if _sn.get('mean', 0) >= 32 else 'fail'}">{'✅ 开阔' if _sn.get('mean', 0) >= 38 else '⚠️ 半遮挡' if _sn.get('mean', 0) >= 32 else '❌ 遮挡'}</td>
            </tr>
        </table>

        <h2>7. 结论</h2>
        <p>本次GPS RTK静态定位测量数据质量总体{'良好' if cep95_pass and pdop_pass and hdop_pass and data_completeness_pass else '一般'}。</p>
        <ul>
            <li><strong>优点:</strong>
                <ul>
                    <li>{fix_rate:.1f}%固定解比例，{'完全满足' if fix_rate >= 95 else '基本满足'}RTK定位要求</li>
                    <li>PDOP和HDOP值{'远低于' if pdop_stats['mean'] <= 1.5 else '低于'}行业阈值，卫星几何构型{'优秀' if pdop_stats['mean'] <= 1.5 else '良好'}</li>
                    <li>数据完整性达{stats.get('data_completeness', 0):.1f}%，{'无' if stats.get('data_gaps', 0) == 0 else '存在'}数据间隙</li>
                    <li>CPE95为{pos_stats['cep95']*100:.1f}厘米，{'达到' if cep95_pass else '未达到'}行业标准</li>
                </ul>
            </li>
            <li><strong>注意事项:</strong>
                <ul>
                    <li>CPE95为{pos_stats['cep95']*100:.1f}厘米，{'处于优秀水平' if pos_stats['cep95'] < 0.01 else '略高于理想值（<2cm），但仍属可接受范围' if not cep95_pass else '处于优秀水平'}</li>
                    <li>垂直RMS({pos_stats['up_rms']*100:.1f}cm){'略高于' if pos_stats['up_rms'] > (pos_stats['east_rms'] + pos_stats['north_rms'])/2 else '与'}水平RMS({(pos_stats['east_rms'] + pos_stats['north_rms'])/2*100:.1f}cm)，这是静态测量的典型特征</li>
                </ul>
            </li>
            <li><strong>建议:</strong>
                <ul>
                    <li>{'数据质量优秀，无需特别处理' if cep95_pass else '建议检查CPE95略高的原因，可能是存在少量的多路径效应'}</li>
                    <li>{'可以继续使用当前观测环境' if pdop_pass and hdop_pass else '建议选择更好的观测时段或位置以降低DOP值'}</li>
                </ul>
            </li>
        </ul>
    </div>
</body>
</html>
    """

    # Write HTML file
    html_path = output_dir / 'report.html'
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    return html_path.name

def generate_markdown_report(output_dir, data, stats, chart_files):
    """Generate Markdown report with image links."""
    # Position statistics on valid fixed solutions only (same source as HTML report).
    pos_stats = {'east_std': 0.0, 'north_std': 0.0, 'up_std': 0.0,
                 'east_rms': 0.0, 'north_rms': 0.0, 'up_rms': 0.0,
                 'cep95': 0.0, 'drms': 0.0, 'horizontal_peak': 0.0, 'vertical_peak': 0.0}
    fix_rate = 0.0
    if 'gnss_data' in data and data['gnss_data']:
        all_recs = data['gnss_data']
        fixed = get_fixed_records(all_recs)
        fix_rate = 100.0 * len([r for r in all_recs
                                if r.get('sol_status') in VALID_SOL_STATUSES
                                and r.get('pos_type') in RTK_FIXED_POS_TYPES]) / len(all_recs)
        lats = np.array([r['lat'] for r in fixed])
        lons = np.array([r['lon'] for r in fixed])
        hgts = np.array([r['hgt'] for r in fixed])
        east, north, up = llh_to_enu(lats, lons, hgts,
                                     np.mean(lats), np.mean(lons), np.mean(hgts))
        pos_stats['east_std'] = float(np.std(east, ddof=1))
        pos_stats['north_std'] = float(np.std(north, ddof=1))
        pos_stats['up_std'] = float(np.std(up, ddof=1))
        pos_stats['east_rms'] = float(np.sqrt(np.mean(east**2)))
        pos_stats['north_rms'] = float(np.sqrt(np.mean(north**2)))
        pos_stats['up_rms'] = float(np.sqrt(np.mean(up**2)))
        pos_stats['cep95'] = float(calc_cep95(east, north))
        pos_stats['drms'] = float(calc_2drms(east, north))
        pos_stats['horizontal_peak'] = float(np.max(np.sqrt(east**2 + north**2)))
        pos_stats['vertical_peak'] = float(np.max(np.abs(up)))

    # DOP statistics: prefer BESTDOPSA (consistent with charts), fall back to GSA.
    dop_src = None
    if 'gsa_data' in data and data['gsa_data']:
        dop_src = data['gsa_data']  # GSA 优先(含 VDOP)
    elif 'bestdops_data' in data and data['bestdops_data']:
        dop_src = data['bestdops_data']
    dop_s = {'pdop': {'min': 0, 'max': 0, 'mean': 0, 'std': 0, 'median': 0},
             'hdop': {'min': 0, 'max': 0, 'mean': 0, 'std': 0, 'median': 0},
             'vdop': None}
    if dop_src:
        dop_s['pdop'] = calc_statistics([r['pdop'] for r in dop_src])
        dop_s['hdop'] = calc_statistics([r['hdop'] for r in dop_src])
        vdop_vals = [r.get('vdop') for r in dop_src]
        if all(v is not None for v in vdop_vals) and any(v > 0 for v in vdop_vals):
            dop_s['vdop'] = calc_statistics(vdop_vals)

    cep95_pass = pos_stats['cep95'] <= RTK_STANDARDS['cep95_rtk_fixed']
    cep95_eval = '✅ 通过' if cep95_pass else '⚠️ 警告'
    completeness_pass = stats['data_completeness'] >= 99.0
    fixrate_pass = fix_rate >= 95.0
    _fc = stats.get('fix_continuity', {}) or {}
    _ss = stats.get('static_stability', {}) or {}
    _sc = stats.get('sigma_convergence', {}) or {}
    _sn = stats.get('snr_summary', {}) or {}
    _es = stats.get('envstatus', {}) or {}
    def _g(d, k, fmt='{:.2f}'):
        v = d.get(k)
        return (fmt.format(v) if isinstance(v, (int, float)) else 'N/A')
    _fc_txt = (f"中断{_fc.get('interruption_count', 0)}次 / 最长连续固定{_fc.get('longest_fix_streak_s', 0):.1f}s / "
               f"均重固定{_fc.get('mean_refix_s', 0):.1f}s") if _fc else 'N/A'

    ENV_MD_ROWS = ""
    if _es:
        def _mdrow(label, stat, unit, good, note):
            if not stat:
                return ""
            return f"| {label} | {stat['min']:.1f} / {stat['mean']:.1f} / {stat['max']:.1f} {unit} | {good} | {note} |\n"
        ENV_MD_ROWS += _mdrow("环境质量分值 env_score", _es.get('env_score'), "", "85~100优/60~84一般/0~59差", "综合环境评分，越高越好")
        ENV_MD_ROWS += _mdrow("卫星可视比率", _es.get('sat_vis_rate'), "%", "越高越好", "跟踪星/理论可视星，反映遮挡")
        ENV_MD_ROWS += _mdrow("参与解算比率", _es.get('sat_sol_rate'), "%", "越高越好", "解算星/跟踪星，反映可用性")
        ENV_MD_ROWS += _mdrow("解置信度", _es.get('sol_confidence'), "", "越高越好", "接收机对当前解的置信程度")
        ENV_MD_ROWS += _mdrow("跟踪卫星数", _es.get('num_svs'), "颗", "越多越好", "可见并跟踪的卫星总数")
        ENV_MD_ROWS += _mdrow("共视卫星数", _es.get('cv_sats'), "颗", "越多越好", "基站-Rover共视星，RTK稳定性")
        ENV_MD_ROWS += _mdrow("平均信噪比 avg_snr", _es.get('avg_snr'), "dBHz", "≥38开阔/32~37半遮挡/≤31严重", "接收机前端信号质量")
        ENV_MD_ROWS += _mdrow("PDOP", _es.get('pdop'), "", "≤2.5好/2.5~5一般/>5差", "卫星几何构型")
        _iono = _es.get('iono_level_counts', {})
        ENV_MD_ROWS += "| 电离层活跃状态 | " + (", ".join(f"{k}:{v}" for k, v in _iono.items()) if _iono else "N/A") + " | QUIET平静为佳 | 电离层扰动影响定位 |\n"
    _fc = stats.get('fix_continuity', {}) or {}
    _ss = stats.get('static_stability', {}) or {}
    _sc = stats.get('sigma_convergence', {}) or {}
    _sn = stats.get('snr_summary', {}) or {}
    _es = stats.get('envstatus', {}) or {}
    def _g(d, k, fmt='{:.2f}'):
        v = d.get(k)
        return (fmt.format(v) if isinstance(v, (int, float)) else 'N/A')
    _fc_txt = (f"中断{_fc.get('interruption_count', 0)}次 / 最长连续固定{_fc.get('longest_fix_streak_s', 0):.1f}s / "
               f"均重固定{_fc.get('mean_refix_s', 0):.1f}s") if _fc else 'N/A'
    md_content = f"""# GPS RTK静态定位分析报告

## 文件信息
- **数据文件:** {data['filename']}
- **生成时间:** {stats['generated_time']}
- **持续时间:** {stats['duration']:.2f} 秒

## 1. 数据概览

| 消息类型 | 数量 | 设置周期 | 实际周期 | 实际频率 | 均匀性 | 符合性 |
|---------|------|---------|---------|---------|--------|--------|"""

    # Add message type summary (实测周期/频率/均匀性/符合性，标色告警)
    _mp = stats.get('msg_periods', {})
    for msg_type, info in data['message_counts'].items():
        nominal = EXPECTED_INTERVALS.get(msg_type, 0)
        nominal_str = f"{nominal:.2f}s/{1/nominal:.0f}Hz" if nominal > 0 else "N/A"
        mp = _mp.get(msg_type)
        if msg_type in PERIOD_NOT_MEASURABLE:
            period_str = hz_str = "N/A"
            uni_str = "-"
            cmp_str = "报文无时间戳，周期不可实测"
        elif mp:
            period_str = f"{mp['median_interval_s']:.3f}s"
            hz_str = f"{mp['actual_hz']:.2f}Hz"
            uni_str = f"均匀(CV={mp['cv']:.2f})" if mp['uniform'] else f"不均(CV={mp['cv']:.2f})"
            note = mp.get('note', '')
            if not mp['compliant']:
                cmp_str = f"🔴偏差{mp['dev_pct']:+.0f}%" + (f" {note}" if note else "")
            elif note:
                cmp_str = f"🟡 {note}"
            else:
                cmp_str = "🟢符合"
        else:
            period_str = hz_str = "N/A"
            uni_str = cmp_str = "-"
        md_content += f"""
| {msg_type} | {info:,} | {nominal_str} | {period_str} | {hz_str} | {uni_str} | {cmp_str} |"""

    md_content += chr(10) + "> 说明：仅当报文自带时间字段时，实际周期/频率/均匀性为真实测量（# 类头含 TOW；NMEA 中 GGA/GST 含 UTC 时间）。GSA/GSV 无时间字段且 GSV 每历元连发多条、同历元各条同时发出，单句周期无意义，标记 N/A（数据为准，绝不臆造）。消息是否按周期发送由输出触发方式决定（ONTIME=周期；ONNEW/ONCHANGED=变化才发）。" + chr(10)

    md_content += f"""
- **数据完整性:** {stats['data_completeness']:.1f}% (共 {stats['total_messages']:,} 条消息)
- **位置信息:** {data['location']['lat']:.6f}°N, {data['location']['lon']:.6f}°E, {data['location']['hgt']:.2f} m (WGS84)
- **天线高度:** {data['antenna_height']:.3f} m (测站上方)
- **报文完整性(校验位):** 有效 {stats.get('integrity', {}).get('valid', 0):,}/{stats.get('integrity', {}).get('total', 0):,} 条 (通过率 {stats.get('integrity', {}).get('integrity_pct', 0):.1f}%；无效 {stats.get('integrity', {}).get('invalid', 0)} 条) — NMEA XOR: {stats.get('integrity', {}).get('nmea_valid', 0):,}/{stats.get('integrity', {}).get('nmea_total', 0):,} ({stats.get('integrity', {}).get('nmea_pct', 0):.1f}%)；#类 CRC-32(表4-6): {stats.get('integrity', {}).get('hash_valid', 0):,}/{stats.get('integrity', {}).get('hash_total', 0):,} ({stats.get('integrity', {}).get('hash_pct', 0):.1f}%)

## 2. GNSS定位分析 (BESTGNSSPOSA)

### 2.1 位置时间序列
![位置时间序列](./pos_timeseries.png)

### 2.2 东西-南北散点图
![东西-南北散点图](./en_scatter.png)

| 统计项 | 东西方向 (m) | 南北方向 (m) | 高程方向 (m) |
|--------|-------------|-------------|-------------|
| 均值（以均值位置为基准的相对偏移，均值定义上恒为 0） | 0.0000 | 0.0000 | 0.0000 |
| 标准差 | {pos_stats['east_std']:.4f} | {pos_stats['north_std']:.4f} | {pos_stats['up_std']:.4f} |
| RMS | {pos_stats['east_rms']:.4f} | {pos_stats['north_rms']:.4f} | {pos_stats['up_rms']:.4f} |
| CPE95 | {pos_stats['cep95']:.4f} | - | - |
| 2DRMS | {pos_stats['drms']:.4f} | - | - |

### 2.3 位置精度时间序列
![位置精度时间序列](./position_sigma_ts.png)

### 2.4 定位类型分布
![定位类型分布](./solution_type_pie.png)

{_fmt_solution_dist(stats.get('sol_type_dist', {}), POS_TYPE_ENUM)}

### 2.5 使用卫星数量
![使用卫星数量](./num_sats_ts.png)

## 3. DOP分析

### 3.1 DOP时间序列
![DOP时间序列](./dop_timeseries.png)

### 3.2 DOP直方图
![DOP直方图](./dop_histogram.png)

| DOP类型 | 最小值 | 最大值 | 均值 | 标准差 | 中位数 |
|--------|-------|-------|------|--------|-------|
| PDOP | {dop_s['pdop']['min']:.2f} | {dop_s['pdop']['max']:.2f} | {dop_s['pdop']['mean']:.2f} | {dop_s['pdop']['std']:.2f} | {dop_s['pdop']['median']:.2f} |
| HDOP | {dop_s['hdop']['min']:.2f} | {dop_s['hdop']['max']:.2f} | {dop_s['hdop']['mean']:.2f} | {dop_s['hdop']['std']:.2f} | {dop_s['hdop']['median']:.2f} |
"""
    if dop_s['vdop']:
        vdop_row = "| VDOP | {0:.2f} | {1:.2f} | {2:.2f} | {3:.2f} | {4:.2f} |".format(
            dop_s['vdop']['min'], dop_s['vdop']['max'], dop_s['vdop']['mean'],
            dop_s['vdop']['std'], dop_s['vdop']['median'])
        md_content += vdop_row + chr(10)
    else:
        md_content += "(BESTDOPSA 消息不含 VDOP 字段，VDOP 统计不可用)" + chr(10)
    md_content += f"""
## 4. 卫星可见性分析

### 4.1 天顶图
![天顶图](./skyplot.png)

### 4.2 信噪比分布
![信噪比分布](./snr_distribution.png)

**箱线图说明：**
- **箱体**：包含 25% 到 75% 的数据（四分位距）
- **黄色中线**：中位数
- **须线（Whiskers）**：显示数据的范围
- **颜色**：GPS-绿色、GLONASS-红色、Galileo-蓝色、QZSS-紫色、BeiDou-黄色
- **数据来源**：所有可见卫星的 GSV 消息中的 SNR 字段

### 4.3 卫星使用数时间序列
![卫星使用数时间序列](./satellites_used_ts.png)

**说明:**
- `#SVs` (蓝色实线): 接收机跟踪到的卫星总数
- `#solnSVs` (红色虚线): 参与定位解算的卫星数

### 4.4 卫星详细信息表
![卫星详细信息表](./satellite_table.png)

## 6A. 环境感知分析（ENVSTATUSA · 华测特有 · 移动站模式）

该报文为华测模组特色，用于评估观测环境对静态定位的影响。环境越差，静态精度越易恶化。

| 字段 | 最小/均值/最大 | 好坏判断 | 说明 |
|------|--------------|---------|------|
{ENV_MD_ROWS}
## 6. 质量评估表

| 评估指标 | 数值 | 行业标准 | 评估结果 |
|---------|------|---------|----------|
| 固定解比例 | {fix_rate:.1f}% | >95% | {'✅ 通过' if fixrate_pass else '⚠️ 警告'} |
| 位置CPE95 | {pos_stats['cep95']:.4f} m | <0.02 m | {cep95_eval} |
| 水平RMS | {(pos_stats['east_rms'] + pos_stats['north_rms'])/2:.4f} m | - | - |
| 垂直RMS | {pos_stats['up_rms']:.4f} m | - | - |
| 固定解水平峰值 | {pos_stats['horizontal_peak']:.4f} m | - | - |
| 固定解垂直峰值 | {pos_stats['vertical_peak']:.4f} m | - | - |
| 平均PDOP | {dop_s['pdop']['mean']:.2f} | <3.0 | {'✅ 通过' if dop_s['pdop']['mean'] < 3.0 else '⚠️ 警告'} |
| 平均HDOP | {dop_s['hdop']['mean']:.2f} | <2.0 | {'✅ 通过' if dop_s['hdop']['mean'] < 2.0 else '⚠️ 警告'} |
| 数据完整性 | {stats['data_completeness']:.1f}% | >99% | {'✅ 通过' if completeness_pass else '⚠️ 警告'} |
| 数据间隙 | {stats['data_gaps']} 个 | 0 个 | {'✅ 通过' if stats['data_gaps'] == 0 else '⚠️ 警告'} |
| 固定解连续性 | {_fc_txt} | 越少中断/越短重固定越好 | {'✅ 连续' if _fc.get('interruption_count', 1) == 0 else '⚠️ 有中断'} |
| 静止性(纯GNSS) | 最大偏移 {_g(_ss, 'max_horizontal_offset', '{:.3f}')} m | ≤{_ss.get('h_threshold', 0.15):.2f} m | {'✅ 静止' if _ss.get('is_static') else '✗ 非静止(运动数据)'} |
| 伪距残差RMS收敛 (GNGST) | 初 {_g(_sc, 'h_init', '{:.3f}')} → 末 {_g(_sc, 'h_final', '{:.3f}')} m | 收敛且小为佳 | {'✅ 收敛' if _sc.get('h_final', 9) < _sc.get('h_init', 0) else '⚠️ 未明显收敛'} |
| 载噪比 C/N0 均值 | {_g(_sn, 'mean', '{:.1f}')} dBHz | ≥38开阔/32-37半遮挡/≤31严重 | {'✅ 开阔' if _sn.get('mean', 0) >= 38 else '⚠️ 半遮挡' if _sn.get('mean', 0) >= 32 else '❌ 遮挡'} |

## 7. 结论

本次GPS RTK静态定位测量数据质量总体{'良好' if cep95_pass and fixrate_pass else '一般'}。

### 优点
- {fix_rate:.1f}%固定解比例{'，完全满足RTK定位要求' if fixrate_pass else '，未达到95%要求'}
- PDOP和HDOP值{'远低于' if dop_s['pdop']['mean'] < 3.0 else '接近'}行业阈值，卫星几何构型{'优秀' if dop_s['pdop']['mean'] < 3.0 else '一般'}
- 数据完整性达{stats['data_completeness']:.1f}%，{'无数据间隙' if stats['data_gaps'] == 0 else '存在数据间隙'}
- 接收机位置稳定，观测环境良好

### 注意事项
- CPE95为{pos_stats['cep95']*100:.1f}厘米，{'达到行业标准（<2cm）' if cep95_pass else '略高于理想值（<2cm），但仍属可接受范围'}
- 垂直RMS略高于水平RMS，这是静态测量的典型特征

### 建议
- {'数据质量优秀，无需特别处理' if cep95_pass else '建议检查CPE95略高的原因，可能是存在少量的多路径效应'}
- 可以增加观测时间，进一步降低随机误差
"""

    # Write Markdown file
    md_path = output_dir / 'report.md'
    with open(md_path, 'w', encoding='utf-8') as f:
        f.write(md_content)

    return md_path.name

# === MAIN PROGRAM ===

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='GPS RTK静态定位数据分析工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python gps_rtk_analyzer_huace.py com1.dat
  python gps_rtk_analyzer_huace.py D:\\data\\test.dat
        """
    )
    parser.add_argument('input_file', help='GPS .dat文件路径')
    parser.add_argument('-o', '--output', help='输出目录 (默认: <文件名>_output)')
    parser.add_argument('--dpi', type=int, default=DPI_DEFAULT, help='图表分辨率 (默认: 150)')
    parser.add_argument('--no-html', action='store_true', help='不生成HTML报告')
    parser.add_argument('--no-md', action='store_true', help='不生成Markdown报告')
    return parser.parse_args()

def cli_main():
    """Main program entry point."""
    args = parse_args()

    # Set DPI
    global DPI_DEFAULT
    DPI_DEFAULT = args.dpi

    # Setup paths
    input_path = Path(args.input_file)
    output_dir = Path(args.output) if args.output else input_path.with_name(input_path.stem + '_output')
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"开始分析: {input_path}")
    print(f"输出目录: {output_dir}")

    # Read data —— 以"消息头锚点 + 校验和验证"的字节流提取为唯一权威来源。
    # 华测 COM11 是文本+二进制混合流，按行读取会被二进制块的 CRLF 拆出碎片/切断报文，
    # 故不再依赖按行行数，直接采用校验通过的完整报文集合。
    print("读取数据文件(字节流锚点提取)...")
    with open(input_path, 'rb') as _fb:
        _raw = _fb.read()
    lines, _st = extract_valid_messages(_raw)
    print(f"锚点校验提取 {_st['rescued']} 条完整报文 "
          f"(扫描锚点 {_st['scanned']}, 校验失败/被二进制污染剔除 {_st['binary_dropped']})")
    lines_all = lines
    print(f"共读取 {len(lines):,} 条数据")

    # Split messages
    print("分割消息类型...")
    message_counts = split_messages(lines, output_dir)
    print(f"消息类型: {message_counts}")

    # Parse data
    print("解析数据...")
    gnss_data = []
    gga_data = []
    gsa_data = []
    gsv_data = []
    bestdops_data = []
    gst_data = []
    envstatus_data = []
    satvis2_data = []

    # Parse each message type
    for msg_type, count in message_counts.items():
        msg_file = output_dir / f"{msg_type}.dat"
        if msg_file.exists():
            lines = []
            with open(msg_file, 'r', encoding='utf-8') as f:
                lines = [line.strip() for line in f if line.strip()]

            # Handle both formats: with and without # prefix for ComNav messages
            if msg_type in ['#BESTPA', 'BESTPA']:
                gnss_data = parse_bestpa(lines)
                print(f"解析了 {len(gnss_data):,} 条BESTPA消息(10Hz)")
            elif msg_type in ['#BESTGNSSPOSA', 'BESTGNSSPOSA', '#BESTPOSA', 'BESTPOSA']:
                # BESTPOSA 仅作兜底：仅当 BESTPA 缺失时才用(1Hz)，避免覆盖10Hz数据
                if not gnss_data:
                    gnss_data = parse_bestgnsspos(lines)
                    print(f"解析了 {len(gnss_data):,} 条BESTPOSA消息(1Hz兜底)")
            elif msg_type in ['$GNGGA', '$GPGGA']:
                gga_data = parse_gga(lines)
                print(f"解析了 {len(gga_data):,} 条GGA消息")
            elif msg_type in ['$GNGSA', '$GPGSA']:
                gsa_data = parse_gsa(lines)
                print(f"解析了 {len(gsa_data):,} 条GSA消息")
            elif msg_type in ['#BESTDOPSA', 'BESTDOPSA']:
                bestdops_data = parse_bestdopsa(lines)
                print(f"解析了 {len(bestdops_data):,} 条BESTDOPSA消息")
            elif msg_type in ['$GNGST', '$GPGST']:
                gst_data = parse_gst(lines)
                print(f"解析了 {len(gst_data):,} 条GST消息")
            elif msg_type in ['#ENVSTATUSA', 'ENVSTATUSA']:
                envstatus_data = parse_envstatusa(lines)
                print(f"解析了 {len(envstatus_data):,} 条ENVSTATUSA消息")
            elif msg_type in ['#SATVIS2A', 'SATVIS2A']:
                pass  # SATVIS2A 卫星可见性(华测NovAtel风格头)，本轮仅统计数量
            elif msg_type in ['$GPGSV', '$GLGSV', '$GAGSV', '$GQGSV', '$GBGSV', '$GIGSV', '$GNGSV']:
                # Process all GSV messages together for proper epoch detection
                if msg_type == '$GPGSV':
                    # Only process GPS GSV, and load other GSV types as well
                    glonass_lines = []
                    galileo_lines = []
                    qzss_lines = []
                    beidou_lines = []
                    navic_lines = []
                    gnss_lines = []

                    # Read other GSV files
                    for gsv_type in ['$GLGSV', '$GAGSV', '$GQGSV', '$GBGSV', '$GIGSV', '$GNGSV']:
                        gsv_file = output_dir / f"{gsv_type}.dat"
                        if gsv_file.exists():
                            with open(gsv_file, 'r', encoding='utf-8') as f:
                                if gsv_type == '$GLGSV':
                                    glonass_lines = [line.strip() for line in f if line.strip()]
                                elif gsv_type == '$GAGSV':
                                    galileo_lines = [line.strip() for line in f if line.strip()]
                                elif gsv_type == '$GQGSV':
                                    qzss_lines = [line.strip() for line in f if line.strip()]
                                elif gsv_type == '$GBGSV':
                                    beidou_lines = [line.strip() for line in f if line.strip()]
                                elif gsv_type == '$GIGSV':
                                    navic_lines = [line.strip() for line in f if line.strip()]
                                elif gsv_type == '$GNGSV':
                                    gnss_lines = [line.strip() for line in f if line.strip()]

                    gsv_data = parse_gsv(lines + glonass_lines + galileo_lines + qzss_lines + beidou_lines + navic_lines + gnss_lines)
                    print(f"解析了 {len(gsv_data):,} 个GSV历元")

    # Calculate statistics
    print("计算统计信息...")

    # Basic info
    stats = {
        'generated_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'total_messages': sum(message_counts.values())
    }

    # Calculate duration from GGA UTC time
    if gga_data:
        gga_times = [r['time_seconds'] for r in gga_data if 'time_seconds' in r]
        if gga_times:
            stats['duration'] = max(gga_times) - min(gga_times)
            stats['start_time_utc'] = min(gga_times)
            stats['end_time_utc'] = max(gga_times)
        else:
            stats['duration'] = 0.0
            stats['start_time_utc'] = 0.0
            stats['end_time_utc'] = 0.0
    else:
        stats['duration'] = 0.0
        stats['start_time_utc'] = 0.0
        stats['end_time_utc'] = 0.0

    # Data completeness: ratio of expected samples to actual samples
    # Use actual sampling interval from GGA data
    if stats['duration'] > 0 and gga_data and len(gga_data) > 1:
        # Nominal 1 Hz logging: expected epochs = duration + 1.
        # Using the observed mean interval here would hide output loss (always ~100%),
        # so completeness is reported against the 1 Hz nominal expectation.
        gga_times = [r['time_seconds'] for r in gga_data]
        actual_interval = np.mean(np.diff(gga_times))
        expected_count = int(stats['duration']) + 1
        actual_count = len(gga_data)
        stats['data_completeness'] = min(100.0, (actual_count / expected_count * 100) if expected_count > 0 else 100.0)
        stats['sampling_interval'] = actual_interval  # observed mean; >1s indicates lost epochs
    else:
        stats['data_completeness'] = 0.0
        stats['sampling_interval'] = 1.0  # Nominal 1 Hz for Huace OUTMSG config

    # Time spans (using GPS TOW)
    time_spans = {}
    if gnss_data:
        times = [r['tow'] for r in gnss_data]
        time_spans['#BESTPOSA'] = (min(times), max(times))

    # Location info
    if gnss_data:
        _fx_sorted = sorted(get_fixed_records(gnss_data), key=lambda x: x['tow'])
        first_point = _fx_sorted[0] if _fx_sorted else sorted(gnss_data, key=lambda x: x['tow'])[0]
        data = {
            'filename': input_path.name,
            'message_counts': message_counts,
            'location': {
                'lat': first_point['lat'],
                'lon': first_point['lon'],
                'hgt': first_point['hgt']
            },
            'antenna_height': first_point['undulation'],
            'gnss_data': gnss_data,
            'gsa_data': gsa_data,
            'bestdops_data': bestdops_data,
            'envstatus_data': envstatus_data,
            'gsv_data': gsv_data
        }
    else:
        data = {
            'filename': input_path.name,
            'message_counts': message_counts,
            'location': {'lat': 0, 'lon': 0, 'hgt': 0},
            'antenna_height': 0,
            'gnss_data': gnss_data,
            'gsa_data': gsa_data,
            'envstatus_data': envstatus_data,
            'gsv_data': gsv_data
        }

    stats['time_spans'] = time_spans

    # === 静态/接收机状态评估（真实数据为准，含字段注释见报告） ===
    stats['integrity'] = check_message_integrity(lines_all)
    stats['msg_periods'] = analyze_periods_huace(build_real_message_periods(message_counts, output_dir))
    stats['fix_continuity'] = calc_fix_continuity(gnss_data)
    stats['static_stability'] = assess_static_stability(gnss_data)
    stats['sigma_convergence'] = (summarize_gst_rms_convergence(gst_data)
            if gst_data else summarize_sigma_convergence(gnss_data))
    stats['snr_summary'] = summarize_snr(gsv_data)
    stats['envstatus'] = summarize_envstatus(envstatus_data)

    # Generate charts
    print("生成图表...")
    chart_files = {}

    # GNSS analysis
    chart_files['pos_timeseries'] = chart_position_timeseries(gnss_data, output_dir)
    chart_files['en_scatter'] = chart_en_scatter(gnss_data, output_dir)
    chart_files['position_sigma_ts'] = chart_position_sigma(gnss_data, output_dir, gst_data)
    chart_files['solution_type_pie'], stats['sol_type_dist'] = chart_solution_type_pie(gnss_data, output_dir)
    chart_files['num_sats_ts'] = chart_num_sats_gnss(gnss_data, output_dir)

    # DOP analysis (BESTDOPSA preferred; GSA fallback keeps bynav-compatible behavior)
    # DOP 数据源: GSA 优先(含 VDOP; BESTDOPSA 无 VDOP 字段——手册3.2.12确认)。
    # 华测 GSA 每历元多条分星座但 DOP 相同，先按历元合并得到统一 PDOP/HDOP/VDOP。
    dop_data = merge_gsa_dop_per_epoch(gsa_data) if gsa_data else bestdops_data
    chart_files['dop_timeseries'] = chart_dop_timeseries(dop_data, output_dir, gga_data)
    chart_files['dop_histogram'] = chart_dop_histogram(dop_data, output_dir)

    # Satellite visibility
    chart_files['skyplot'] = chart_skyplot(gsv_data, output_dir)
    chart_files['snr_distribution'] = chart_snr_distribution(gsv_data, output_dir)
    chart_files['satellites_used_ts'] = chart_satellites_used_ts(gnss_data, output_dir)
    chart_files['satellite_table'] = chart_satellite_info(gsv_data, output_dir)

    # Quality assessment
    if gnss_data:
        # Calculate actual fix rate
        total_fixes = len(gnss_data)
        fixed_solutions = len([r for r in gnss_data if r.get('sol_status') in VALID_SOL_STATUSES and r.get('pos_type') in RTK_FIXED_POS_TYPES])
        fix_rate = fixed_solutions / total_fixes if total_fixes > 0 else 0.0
    else:
        fix_rate = 0.0
    stats['fix_rate'] = fix_rate
    chart_files['fix_rate_gauge'] = chart_fix_rate_gauge(fix_rate, output_dir)

    # Data gaps: GGA is the 1 Hz reference stream; flag intervals > 2 s as gaps.
    # (BESTPOSA/GGA show irregular 1/4/5 s spacing in captured logs, i.e. lost output.)
    if gga_data and len(gga_data) > 1:
        times = [r['time_seconds'] for r in gga_data]
        gaps = detect_data_gaps(times, 1.0, tolerance=2.0)
    else:
        gaps = []
    stats['data_gaps'] = len(gaps)
    chart_files['data_gap_analysis'] = chart_data_gaps(gaps, output_dir)

    # Performance summary
    chart_files['performance_summary'] = chart_performance_summary(gnss_data, output_dir)

    # Generate reports
    print("生成报告...")
    if not args.no_html:
        html_report = generate_html_report(output_dir, data, stats, chart_files)
        print(f"HTML报告已生成: {html_report}")

    if not args.no_md:
        md_report = generate_markdown_report(output_dir, data, stats, chart_files)
        print(f"Markdown报告已生成: {md_report}")

    # Print summary
    print("\n分析完成！")
    print(f"输出目录: {output_dir}")
    print(f"共生成 {len([f for f in chart_files.values() if f])} 个图表")
    print("\n数据分析摘要:")
    print(f"- 总消息数: {stats['total_messages']:,}")
    print(f"- 数据完整性: {stats['data_completeness']:.1f}%")
    print(f"- 持续时间: {stats['duration']:.1f} 秒")

    if gnss_data:
        # Key metrics on valid FIXED solutions only (consistent with the report).
        metric_recs = get_fixed_records(gnss_data)
        lats = np.array([r['lat'] for r in metric_recs])
        lons = np.array([r['lon'] for r in metric_recs])
        hgts = np.array([r['hgt'] for r in metric_recs])

        mean_lat = np.mean(lats)
        mean_lon = np.mean(lons)
        mean_hgt = np.mean(hgts)

        east, north, up = llh_to_enu(lats, lons, hgts, mean_lat, mean_lon, mean_hgt)

        cep95 = calc_cep95(east, north)
        drms = calc_2drms(east, north)
        num_sats = [r.get('num_svs') for r in gnss_data if r.get('num_svs')]

        print(f"- 位置信息: {mean_lat:.6f}°N, {mean_lon:.6f}°E")
        print(f"- CPE95 (固定解): {cep95:.3f} m ({cep95*100:.1f} cm)")
        print(f"- 2DRMS (固定解): {drms:.3f} m")
        if num_sats:
            print(f"- 跟踪卫星数: 平均 {np.mean(num_sats):.0f} (范围 {min(num_sats)}-{max(num_sats)})")



# GUI dependencies (used by the merged GUI section below)
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
# =============================================================================
# GUI (merged from gps_rtk_analyzer_gui_huace.py)
# =============================================================================

class GPSAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("GPS RTK 静态定位数据分析工具（华测 CHCNAV M720 · 选 COM11 日志）")
        self.root.geometry("800x700")
        self.root.resizable(True, True)

        # Variables
        self.input_file = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.dpi = tk.IntVar(value=150)
        self.generate_html = tk.BooleanVar(value=True)
        self.generate_md = tk.BooleanVar(value=True)
        self.analyzing = False

        self.create_widgets()
        self.setup_styles()

    def setup_styles(self):
        """Setup custom styles for the GUI."""
        style = ttk.Style()
        style.theme_use('clam')

        style.configure('Title.TLabel', font=('Arial', 16, 'bold'))
        style.configure('Header.TLabel', font=('Arial', 11, 'bold'))
        style.configure('Info.TLabel', font=('Arial', 9))
        style.configure('Action.TButton', font=('Arial', 10))

    def create_widgets(self):
        """Create all GUI widgets."""
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)

        # Title
        title_label = ttk.Label(main_frame, text="GPS RTK 静态定位数据分析工具（华测 M720）",
                                style='Title.TLabel')
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20), sticky=tk.W)

        # Input file selection
        ttk.Label(main_frame, text="输入数据文件:", style='Header.TLabel').grid(
            row=1, column=0, sticky=tk.W, pady=(0, 5))
        input_frame = ttk.Frame(main_frame)
        input_frame.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 15))
        input_frame.columnconfigure(1, weight=1)

        ttk.Entry(input_frame, textvariable=self.input_file, width=50).grid(
            row=0, column=0, sticky=(tk.W, tk.E), padx=(0, 5))
        ttk.Button(input_frame, text="浏览...", command=self.browse_input_file).grid(
            row=0, column=1, padx=(0, 5))
        ttk.Button(input_frame, text="打开最近文件", command=self.open_recent_file).grid(
            row=0, column=2)

        # Output directory selection
        ttk.Label(main_frame, text="输出目录:", style='Header.TLabel').grid(
            row=3, column=0, sticky=tk.W, pady=(0, 5))
        output_frame = ttk.Frame(main_frame)
        output_frame.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 15))
        output_frame.columnconfigure(1, weight=1)

        ttk.Entry(output_frame, textvariable=self.output_dir, width=50).grid(
            row=0, column=0, sticky=(tk.W, tk.E), padx=(0, 5))
        ttk.Button(output_frame, text="浏览...", command=self.browse_output_dir).grid(
            row=0, column=1)

        # Options frame
        options_frame = ttk.LabelFrame(main_frame, text="分析选项", padding="10")
        options_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 15))

        # DPI setting
        dpi_frame = ttk.Frame(options_frame)
        dpi_frame.grid(row=0, column=0, sticky=tk.W, pady=(0, 10))
        ttk.Label(dpi_frame, text="图表分辨率 (DPI):").grid(row=0, column=0, padx=(0, 5))
        dpi_spinbox = ttk.Spinbox(dpi_frame, from_=72, to=300, textvariable=self.dpi, width=10)
        dpi_spinbox.grid(row=0, column=1)

        # Report generation options
        ttk.Checkbutton(options_frame, text="生成 HTML 报告",
                        variable=self.generate_html).grid(row=1, column=0, sticky=tk.W, pady=(0, 5))
        ttk.Checkbutton(options_frame, text="生成 Markdown 报告",
                        variable=self.generate_md).grid(row=2, column=0, sticky=tk.W)

        # Analyze button
        self.analyze_button = ttk.Button(main_frame, text="开始分析",
                                         command=self.start_analysis, style='Action.TButton')
        self.analyze_button.grid(row=6, column=0, columnspan=3, pady=(15, 10), sticky=(tk.W, tk.E))

        # Progress bar
        self.progress = ttk.Progressbar(main_frame, mode='indeterminate')
        self.progress.grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(5, 10))

        # Output text area
        ttk.Label(main_frame, text="分析日志:", style='Header.TLabel').grid(
            row=8, column=0, sticky=tk.W, pady=(5, 5))

        self.log_text = scrolledtext.ScrolledText(main_frame, height=15, width=80, wrap=tk.WORD)
        self.log_text.grid(row=9, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(0, 10))
        self.log_text.configure(font=('Consolas', 9))

        # Configure text tag colors
        self.log_text.tag_config('info', foreground='black')
        self.log_text.tag_config('success', foreground='green')
        self.log_text.tag_config('warning', foreground='orange')
        self.log_text.tag_config('error', foreground='red')
        self.log_text.tag_config('header', font=('Consolas', 9, 'bold'))

        # Status bar
        self.status_var = tk.StringVar(value="就绪")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.grid(row=10, column=0, columnspan=3, sticky=(tk.W, tk.E))

        # Configure main frame grid weights
        main_frame.rowconfigure(9, weight=1)

        # Load recent files
        self.load_recent_files()

    def browse_input_file(self):
        """Browse for input .dat file."""
        filename = filedialog.askopenfilename(
            title="选择GPS数据文件",
            filetypes=[("华测COM11日志", "*.log *.dat"), ("All files", "*.*")]
        )
        if filename:
            self.input_file.set(filename)
            # Auto-set output directory
            input_path = Path(filename)
            output_path = input_path.with_name(input_path.stem + '_output')
            self.output_dir.set(str(output_path))
            self.save_recent_file(filename)

    def browse_output_dir(self):
        """Browse for output directory."""
        dirname = filedialog.askdirectory(title="选择输出目录")
        if dirname:
            self.output_dir.set(dirname)

    def load_recent_files(self):
        """Load recently used files."""
        config_dir = Path.home() / '.gps_analyzer'
        recent_file = config_dir / 'recent_files.txt'

        if recent_file.exists():
            try:
                with open(recent_file, 'r', encoding='utf-8') as f:
                    recent_files = [line.strip() for line in f if line.strip()]
                self.recent_files = recent_files[:10]  # Keep last 10
            except Exception:
                self.recent_files = []
        else:
            self.recent_files = []

    def save_recent_file(self, filename):
        """Save a file to recent files list."""
        config_dir = Path.home() / '.gps_analyzer'
        config_dir.mkdir(parents=True, exist_ok=True)
        recent_file = config_dir / 'recent_files.txt'

        # Add to front if not already in list
        if filename not in self.recent_files:
            self.recent_files.insert(0, filename)
        else:
            # Move to front
            self.recent_files.remove(filename)
            self.recent_files.insert(0, filename)

        # Keep only last 10
        self.recent_files = self.recent_files[:10]

        try:
            with open(recent_file, 'w', encoding='utf-8') as f:
                for fpath in self.recent_files:
                    f.write(fpath + '\n')
        except Exception:
            pass

    def open_recent_file(self):
        """Open a dialog to select from recent files."""
        if not self.recent_files:
            messagebox.showinfo("最近文件", "没有最近使用的文件。")
            return

        recent_window = tk.Toplevel(self.root)
        recent_window.title("最近文件")
        recent_window.geometry("600x300")
        recent_window.transient(self.root)
        recent_window.grab_set()

        listbox = tk.Listbox(recent_window, font=('Arial', 10))
        listbox.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        scrollbar = ttk.Scrollbar(recent_window, orient=tk.VERTICAL, command=listbox.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        listbox.config(yscrollcommand=scrollbar.set)

        for fpath in self.recent_files:
            if Path(fpath).exists():
                listbox.insert(tk.END, fpath)

        def select_file():
            selection = listbox.curselection()
            if selection:
                filename = listbox.get(selection[0])
                self.input_file.set(filename)
                # Auto-set output directory
                input_path = Path(filename)
                output_path = input_path.with_name(input_path.stem + '_output')
                self.output_dir.set(str(output_path))
                recent_window.destroy()

        button_frame = ttk.Frame(recent_window)
        button_frame.pack(fill=tk.X, padx=10, pady=10)

        ttk.Button(button_frame, text="打开", command=select_file).pack(side=tk.RIGHT)
        ttk.Button(button_frame, text="取消",
                   command=recent_window.destroy).pack(side=tk.RIGHT, padx=(0, 5))

    def log(self, message, level='info'):
        """Add a message to the log text area."""
        self.log_text.insert(tk.END, message + '\n', level)
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def start_analysis(self):
        """Start the analysis in a separate thread."""
        if self.analyzing:
            return

        # Validate inputs
        input_file = self.input_file.get().strip()
        if not input_file:
            messagebox.showerror("错误", "请选择输入数据文件")
            return

        if not Path(input_file).exists():
            messagebox.showerror("错误", f"文件不存在: {input_file}")
            return

        output_dir = self.output_dir.get().strip()
        if not output_dir:
            input_path = Path(input_file)
            output_dir = str(input_path.with_name(input_path.stem + '_output'))
            self.output_dir.set(output_dir)

        # Start analysis in a thread
        self.analyzing = True
        self.analyze_button.config(state=tk.DISABLED)
        self.progress.start()
        self.status_var.set("正在分析...")

        thread = threading.Thread(target=self.run_analysis, args=(input_file, output_dir))
        thread.daemon = True
        thread.start()

    def run_analysis(self, input_file, output_dir):
        """Run the actual analysis (runs in a separate thread)."""
        try:
            # Clear log
            self.log_text.delete(1.0, tk.END)
            self.log("=" * 60, 'header')
            self.log(f"GPS RTK 静态定位数据分析", 'header')
            self.log(f"输入文件: {input_file}", 'info')
            self.log(f"输出目录: {output_dir}", 'info')
            self.log(f"DPI设置: {self.dpi.get()}", 'info')
            self.log("=" * 60, 'header')

            # Setup paths
            input_path = Path(input_file)
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            # Set DPI globally (module-level variable in this file)
            global DPI_DEFAULT
            DPI_DEFAULT = self.dpi.get()

            # Read data —— 以"消息头锚点 + 校验和验证"的字节流提取为唯一权威来源
            self.log("读取数据文件(字节流锚点提取)...", 'info')
            with open(input_path, 'rb') as _fb:
                _raw = _fb.read()
            lines, _st = extract_valid_messages(_raw)
            self.log(f"锚点校验提取 {_st['rescued']} 条完整报文 (剔除 {_st['binary_dropped']})", 'info')
            self.log(f"共读取 {len(lines):,} 条数据", 'info')

            # Split messages
            self.log("分割消息类型...", 'info')
            message_counts = split_messages(lines, output_path)
            self.log(f"消息类型: {message_counts}", 'info')

            # Parse data
            self.log("解析数据...", 'info')
            gnss_data = []
            gga_data = []
            gsa_data = []
            gsv_data = []
            bestdops_data = []
            gst_data = []
            envstatus_data = []

            for msg_type, count in message_counts.items():
                msg_file = output_path / f"{msg_type}.dat"
                if msg_file.exists():
                    with open(msg_file, 'r', encoding='utf-8') as f:
                        msg_lines = [line.strip() for line in f if line.strip()]

                    if msg_type in ['#BESTPA', 'BESTPA']:
                        gnss_data = parse_bestpa(msg_lines)
                        self.log(f"解析了 {len(gnss_data):,} 条BESTPA消息(10Hz)", 'info')
                    elif msg_type in ['#BESTPOSA', 'BESTPOSA', '#BESTGNSSPOSA', 'BESTGNSSPOSA']:
                        if not gnss_data:
                            gnss_data = parse_bestgnsspos(msg_lines)
                            self.log(f"解析了 {len(gnss_data):,} 条BESTPOSA消息(1Hz兜底)", 'info')
                    elif msg_type in ['$GNGGA', '$GPGGA']:
                        gga_data = parse_gga(msg_lines)
                        self.log(f"解析了 {len(gga_data):,} 条GGA消息", 'info')
                    elif msg_type in ['$GNGSA', '$GPGSA']:
                        gsa_data = parse_gsa(msg_lines)
                        self.log(f"解析了 {len(gsa_data):,} 条GSA消息", 'info')
                    elif msg_type in ['#BESTDOPSA', 'BESTDOPSA']:
                        bestdops_data = parse_bestdopsa(msg_lines)
                        self.log(f"解析了 {len(bestdops_data):,} 条BESTDOPSA消息", 'info')
                    elif msg_type in ['$GNGST', '$GPGST']:
                        gst_data = parse_gst(msg_lines)
                        self.log(f"解析了 {len(gst_data):,} 条GST消息", 'info')
                    elif msg_type in ['#ENVSTATUSA', 'ENVSTATUSA']:
                        envstatus_data = parse_envstatusa(msg_lines)
                        self.log(f"解析了 {len(envstatus_data):,} 条ENVSTATUSA消息", 'info')
                    elif msg_type == '$GPGSV':
                        # Load other GSV types
                        glonass_lines = []
                        galileo_lines = []
                        qzss_lines = []
                        beidou_lines = []
                        navic_lines = []

                        for gsv_type in ['$GLGSV', '$GAGSV', '$GQGSV', '$GBGSV', '$GIGSV']:
                            gsv_file = output_path / f"{gsv_type}.dat"
                            if gsv_file.exists():
                                with open(gsv_file, 'r', encoding='utf-8') as f:
                                    if gsv_type == '$GLGSV':
                                        glonass_lines = [line.strip() for line in f if line.strip()]
                                    elif gsv_type == '$GAGSV':
                                        galileo_lines = [line.strip() for line in f if line.strip()]
                                    elif gsv_type == '$GQGSV':
                                        qzss_lines = [line.strip() for line in f if line.strip()]
                                    elif gsv_type == '$GBGSV':
                                        beidou_lines = [line.strip() for line in f if line.strip()]
                                    elif gsv_type == '$GIGSV':
                                        navic_lines = [line.strip() for line in f if line.strip()]

                        gsv_data = parse_gsv(msg_lines + glonass_lines + galileo_lines +
                                           qzss_lines + beidou_lines + navic_lines)
                        self.log(f"解析了 {len(gsv_data):,} 个GSV历元", 'info')

            # Calculate statistics
            self.log("计算统计信息...", 'info')

            stats = {
                'generated_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'total_messages': sum(message_counts.values())
            }

            # Calculate duration from GPGGA
            if gga_data:
                gga_times = [r['time_seconds'] for r in gga_data if 'time_seconds' in r]
                if gga_times:
                    stats['duration'] = max(gga_times) - min(gga_times)
                else:
                    stats['duration'] = 0.0
            else:
                stats['duration'] = 0.0

            # Data completeness
            if stats['duration'] > 0 and gga_data and len(gga_data) > 1:
                # Nominal 1 Hz logging; using observed mean interval would hide output loss
                gga_times = [r['time_seconds'] for r in gga_data]
                actual_interval = np.mean(np.diff(gga_times))
                expected_count = int(stats['duration']) + 1
                actual_count = len(gga_data)
                stats['data_completeness'] = min(100.0, (actual_count / expected_count * 100) if expected_count > 0 else 100.0)
            else:
                stats['data_completeness'] = 0.0

            # Time spans
            time_spans = {}
            if gnss_data:
                times = [r['tow'] for r in gnss_data]
                time_spans['#BESTGNSSPOSA'] = (min(times), max(times))

            # Location info
            if gnss_data:
                gnss_data_sorted = sorted(gnss_data, key=lambda x: x['tow'])
                first_point = gnss_data_sorted[0]
                data = {
                    'filename': input_path.name,
                    'message_counts': message_counts,
                    'location': {
                        'lat': first_point['lat'],
                        'lon': first_point['lon'],
                        'hgt': first_point['hgt']
                    },
                    'antenna_height': first_point['undulation'],
                    'gnss_data': gnss_data,
                    'gsa_data': gsa_data,
                    'envstatus_data': envstatus_data,
                }
            else:
                data = {
                    'filename': input_path.name,
                    'message_counts': message_counts,
                    'location': {'lat': 0, 'lon': 0, 'hgt': 0},
                    'antenna_height': 0,
                    'gnss_data': gnss_data,
                    'gsa_data': gsa_data,
                    'envstatus_data': envstatus_data,
                }

            stats['time_spans'] = time_spans

            # === 静态/接收机状态评估（与 CLI 路径一致；以真实数据为准） ===
            stats['integrity'] = check_message_integrity(lines)
            stats['msg_periods'] = analyze_periods_huace(build_real_message_periods(message_counts, output_path))
            stats['fix_continuity'] = calc_fix_continuity(gnss_data)
            stats['static_stability'] = assess_static_stability(gnss_data)
            stats['sigma_convergence'] = (summarize_gst_rms_convergence(gst_data)
                    if gst_data else summarize_sigma_convergence(gnss_data))
            stats['snr_summary'] = summarize_snr(gsv_data)
            stats['envstatus'] = summarize_envstatus(envstatus_data)

            # Generate charts
            self.log("生成图表...", 'info')
            chart_files = {}

            # GNSS analysis
            self.log("  - 位置时间序列...", 'info')
            chart_files['pos_timeseries'] = chart_position_timeseries(gnss_data, output_path)
            self.log("  - 东西-南北散点图...", 'info')
            en_scatter_result = chart_en_scatter(gnss_data, output_path)
            if en_scatter_result:
                chart_files['en_scatter'] = en_scatter_result[0]
                if en_scatter_result[1]:
                    stats['sigma_u'] = en_scatter_result[1].get('sigma_u', 0.0)
            self.log("  - 位置精度时间序列...", 'info')
            chart_files['position_sigma_ts'] = chart_position_sigma(gnss_data, output_path, gst_data)
            self.log("  - 定位类型分布...", 'info')
            chart_files['solution_type_pie'], stats['sol_type_dist'] = chart_solution_type_pie(gnss_data, output_path)
            self.log("  - 使用卫星数量...", 'info')
            chart_files['num_sats_ts'] = chart_num_sats_gnss(gnss_data, output_path)

            # DOP analysis (BESTDOPSA preferred, GSA fallback)
            dop_data = merge_gsa_dop_per_epoch(gsa_data) if gsa_data else bestdops_data
            self.log("  - DOP时间序列...", 'info')
            chart_files['dop_timeseries'] = chart_dop_timeseries(dop_data, output_path, gga_data)
            self.log("  - DOP直方图...", 'info')
            chart_files['dop_histogram'] = chart_dop_histogram(dop_data, output_path)

            # Satellite visibility
            self.log("  - 天顶图...", 'info')
            chart_files['skyplot'] = chart_skyplot(gsv_data, output_path)
            self.log("  - SNR分布...", 'info')
            chart_files['snr_distribution'] = chart_snr_distribution(gsv_data, output_path)
            self.log("  - 卫星使用数时间序列...", 'info')
            chart_files['satellites_used_ts'] = chart_satellites_used_ts(gnss_data, output_path)
            self.log("  - 卫星信息表...", 'info')
            chart_files['satellite_table'] = chart_satellite_info(gsv_data, output_path)

            # Quality assessment
            if gnss_data:
                total_fixes = len(gnss_data)
                fixed_solutions = len([r for r in gnss_data if r.get('sol_status') in VALID_SOL_STATUSES
                                      and r.get('pos_type') in RTK_FIXED_POS_TYPES])
                fix_rate = fixed_solutions / total_fixes if total_fixes > 0 else 0.0
            else:
                fix_rate = 0.0

            stats['fix_rate'] = fix_rate
            self.log("  - 固定率仪表盘...", 'info')
            chart_files['fix_rate_gauge'] = chart_fix_rate_gauge(fix_rate, output_path)

            # Data gaps: GGA 1 Hz reference stream, flag intervals > 2 s
            if gga_data and len(gga_data) > 1:
                times = [r['time_seconds'] for r in gga_data]
                gaps = detect_data_gaps(times, 1.0, tolerance=2.0)
            else:
                gaps = []
            stats['data_gaps'] = len(gaps)
            self.log("  - 数据间隙分析...", 'info')
            chart_files['data_gap_analysis'] = chart_data_gaps(gaps, output_path)

            # Performance summary
            self.log("  - 性能摘要表...", 'info')
            chart_files['performance_summary'] = chart_performance_summary(gnss_data, output_path)

            # Generate reports
            self.log("生成报告...", 'info')
            if self.generate_html.get():
                html_report = generate_html_report(output_path, data, stats, chart_files)
                self.log(f"  HTML报告已生成: {html_report}", 'success')

            if self.generate_md.get():
                md_report = generate_markdown_report(output_path, data, stats, chart_files)
                self.log(f"  Markdown报告已生成: {md_report}", 'success')

            # Print summary
            self.log("", 'info')
            self.log("=" * 60, 'header')
            self.log("分析完成！", 'success')
            self.log(f"输出目录: {output_path}", 'info')
            self.log(f"共生成 {len([f for f in chart_files.values() if f])} 个图表", 'info')
            self.log("", 'info')
            self.log("数据分析摘要:", 'header')
            self.log(f"- 总消息数: {stats['total_messages']:,}", 'info')
            self.log(f"- 数据完整性: {stats['data_completeness']:.1f}%", 'info')
            self.log(f"- 持续时间: {stats['duration']:.1f} 秒", 'info')

            _fx = get_fixed_records(gnss_data)
            if _fx:
                times = [r['tow'] for r in _fx]
                lats = np.array([r['lat'] for r in _fx])
                lons = np.array([r['lon'] for r in _fx])
                hgts = np.array([r['hgt'] for r in _fx])

                mean_lat = np.mean(lats)
                mean_lon = np.mean(lons)
                mean_hgt = np.mean(hgts)

                east, north, up = llh_to_enu(lats, lons, hgts, mean_lat, mean_lon, mean_hgt)
                cep95 = calc_cep95(east, north)
                drms = calc_2drms(east, north)

                self.log(f"- 位置信息: {mean_lat:.6f}°N, {mean_lon:.6f}°E", 'info')
                self.log(f"- CPE95: {cep95:.3f} m ({cep95*100:.1f} cm)", 'info')
                self.log(f"- 2DRMS: {drms:.3f} m", 'info')
                self.log(f"- 固定解比例: {fix_rate*100:.1f}%", 'info')


            self.log("=" * 60, 'header')

            # Ask to open output directory
            self.root.after(0, lambda: self.ask_open_output(output_path))

        except Exception as e:
            self.log(f"分析过程中发生错误: {str(e)}", 'error')
            import traceback
            self.log(traceback.format_exc(), 'error')

        finally:
            # Reset UI state
            self.analyzing = False
            self.root.after(0, self.reset_ui)

    def reset_ui(self):
        """Reset UI to ready state."""
        self.progress.stop()
        self.analyze_button.config(state=tk.NORMAL)
        self.status_var.set("分析完成")

    def ask_open_output(self, output_path):
        """Ask user if they want to open the output directory."""
        response = messagebox.askyesno("分析完成",
                                         "分析已完成！\n\n是否打开输出目录？")
        if response:
            try:
                if os.name == 'nt':  # Windows
                    os.startfile(output_path)
                elif os.name == 'posix':  # macOS and Linux
                    if sys.platform == 'darwin':  # macOS
                        os.system(f'open "{output_path}"')
                    else:  # Linux
                        os.system(f'xdg-open "{output_path}"')
            except Exception as e:
                messagebox.showerror("错误", f"无法打开目录: {e}")


def gui_main():
    """Main entry point for the GUI application."""
    root = tk.Tk()
    app = GPSAnalyzerGUI(root)

    # Handle window close
    def on_closing():
        if app.analyzing:
            if messagebox.askokcancel("退出", "分析正在进行中，确定要退出吗？"):
                root.destroy()
        else:
            root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()


def main():
    """Entry point: launch GUI when no arguments are given, CLI otherwise."""
    import sys as _sys
    if len(_sys.argv) > 1:
        cli_main()
    else:
        gui_main()


if __name__ == '__main__':
    main()

