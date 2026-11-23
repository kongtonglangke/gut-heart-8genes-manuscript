#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Chunked HTTP-Range downloader.

GCS blocks plain GET (0 bytes, exit 35) on this network but honours Range
requests, so we pull the file in small byte slices and reassemble.

Usage:
  python 11_range_download.py <remote_url> <local_path> [chunk_bytes]
"""
import os
import sys
import time
import urllib.request

UA = {"User-Agent": "Mozilla/5.0"}


def total_size(url):
    req = urllib.request.Request(url, headers=UA, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        cl = r.headers.get("Content-Length")
    return int(cl) if cl else None


def fetch(url, start, end, tries=8):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, headers=dict(UA, Range="bytes=%d-%d" % (start, end)))
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:
            last = e
            time.sleep(1.5 * (i + 1))
    raise last


def main():
    url = sys.argv[1]
    out = sys.argv[2]
    chunk = int(sys.argv[3]) if len(sys.argv) > 3 else 65536

    if os.path.exists(out) and os.path.getsize(out) > 0:
        print("exists, skip:", out, os.path.getsize(out))
        return 0

    n = total_size(url)
    if not n:
        print("cannot determine size")
        return 1
    print("remote size = %d bytes, chunk = %d -> %d requests" % (n, chunk, (n + chunk - 1) // chunk))

    tmp = out + ".part"
    t0 = time.time()
    done = 0
    with open(tmp, "wb") as f:
        while done < n:
            a = done
            b = min(done + chunk - 1, n - 1)
            buf = fetch(url, a, b)
            f.write(buf)
            done += len(buf)
            if done % (chunk * 5) < chunk or done == n:
                el = time.time() - t0
                print("  %6.1f%%  %d/%d  %.0fs  %.1f KB/s" % (
                    100.0 * done / n, done, n, el, done / 1024.0 / max(el, 1e-6)))
            if len(buf) == 0:
                print("empty chunk, abort")
                return 1
    os.replace(tmp, out)
    print("done %s  %d bytes  %.1fs" % (out, n, time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
