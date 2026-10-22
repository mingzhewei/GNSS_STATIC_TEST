"""Focused regressions for both checked-in Huace analyzer copies."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from typing import Optional


ROOT = Path(__file__).resolve().parents[1]


def _load_analyzer(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load analyzer at {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ANALYZERS = {
    "root": _load_analyzer("huace_analyzer_root", ROOT / "gps_rtk_analyzer_huace.py"),
    "Build": _load_analyzer(
        "huace_analyzer_build", ROOT / "Build" / "gps_rtk_analyzer_huace.py"
    ),
}


def _valid_gga_message(analyzer) -> bytes:
    body = b"GNGGA,123519,4807.038,N,01131.000,E,1,08,0.9,545.4,M,46.9,M,,"
    checksum = analyzer._nmea_checksum(body)
    return b"$" + body + b"*" + f"{checksum:02X}".encode("ascii")


def _valid_bestpa_message(analyzer) -> bytes:
    body = (
        b"BESTPA,46.5,0,0,COM1,FINESTEERING,2437,12092300,ab0acbc1,00000200;"
        b"SOL_COMPUTED,NARROW_INT,61,0,31.2304,121.4737,10.0,0.0,"
        b"0.01,0.02,0.03,0,0,20,15,0,0,0,0,0,0,0,0042"
    )
    checksum = analyzer._crc32_block(body)
    return b"#" + body + b"*" + f"{checksum:08x}".encode("ascii")


def _gga_record(station_id: Optional[str] = None) -> str:
    fields = [
        "$GNGGA", "123519", "4807.038", "N", "01131.000", "E", "1", "08",
        "0.9", "545.4", "M", "46.9", "M", "1.0",
    ]
    if station_id is not None:
        fields.append(station_id)
    return ",".join(fields)


class HuaceAnalyzerRegressionTests(unittest.TestCase):
    def test_known_checksums_recover_after_each_single_byte_prefix(self):
        for copy_name, analyzer in ANALYZERS.items():
            with self.subTest(copy=copy_name):
                messages = (
                    _valid_gga_message(analyzer),
                    _valid_bestpa_message(analyzer),
                )
                for message in messages:
                    for prefix in range(256):
                        with self.subTest(
                            copy=copy_name,
                            head=message[:7].decode("ascii", "replace"),
                            prefix=prefix,
                        ):
                            extracted, stats = analyzer.extract_valid_messages(
                                bytes((prefix,)) + message
                            )
                            self.assertEqual(extracted, [message.decode("ascii")])
                            self.assertEqual(stats["rescued"], 1)

    def test_gga_parser_keeps_optional_station_id_behavior(self):
        for copy_name, analyzer in ANALYZERS.items():
            with self.subTest(copy=copy_name):
                without_station = analyzer.parse_gga([_gga_record()])
                with_station = analyzer.parse_gga([_gga_record("0042")])
                malformed_short = analyzer.parse_gga(["$GNGGA,123519,4807.038"])

                self.assertEqual(len(without_station), 1)
                self.assertEqual(without_station[0]["station_id"], "")
                self.assertEqual(len(with_station), 1)
                self.assertEqual(with_station[0]["station_id"], "0042")
                self.assertEqual(malformed_short, [])

    def test_position_sigma_availability_uses_only_valid_gnss_sigma_rows(self):
        for copy_name, analyzer in ANALYZERS.items():
            with self.subTest(copy=copy_name):
                valid_gst = [{
                    "time_seconds": 1.0,
                    "lat_sigma": 0.1,
                    "lon_sigma": 0.2,
                    "alt_sigma": 0.3,
                    "rms": 0.2,
                }]
                invalid_gnss = [
                    {"tow": 1.0, "lat_sigma": None, "lon_sigma": None, "hgt_sigma": None},
                    {"tow": 2.0, "lat_sigma": "bad", "lon_sigma": 0.2, "hgt_sigma": 0.3},
                    {"lat_sigma": 0.01, "lon_sigma": 0.02, "hgt_sigma": 0.03},
                ]
                self.assertEqual(
                    analyzer._position_sigma_availability(invalid_gnss, valid_gst),
                    "gngst-fallback",
                )

                constant_gnss = [
                    {"tow": 1.0, "lat_sigma": 0.01, "lon_sigma": 0.02, "hgt_sigma": 0.03},
                    {"tow": 2.0, "lat_sigma": 0.01, "lon_sigma": 0.02, "hgt_sigma": 0.03},
                ]
                self.assertEqual(
                    analyzer._position_sigma_availability(
                        invalid_gnss + constant_gnss, valid_gst
                    ),
                    "constant-output",
                )

                varied_gnss = [
                    {"tow": 1.0, "lat_sigma": 0.01, "lon_sigma": 0.02, "hgt_sigma": 0.03},
                    {"tow": 2.0, "lat_sigma": 0.02, "lon_sigma": 0.02, "hgt_sigma": 0.03},
                ]
                self.assertEqual(
                    analyzer._position_sigma_availability(varied_gnss, valid_gst),
                    "parsed",
                )
                self.assertEqual(analyzer._position_sigma_availability([], []), "unavailable")
                self.assertEqual(
                    analyzer._position_sigma_availability(
                        [], [{"time_seconds": 1.0, "rms": 0.2}]
                    ),
                    "unavailable",
                )
                self.assertEqual(
                    analyzer._position_sigma_availability([], valid_gst), "gngst-fallback"
                )
                self.assertEqual(
                    analyzer._position_sigma_availability(
                        [], [dict(valid_gst[0], lat_sigma="bad")]
                    ),
                    "unavailable",
                )


if __name__ == "__main__":
    unittest.main()
