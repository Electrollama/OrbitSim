#!/usr/bin/env python3
"""Test runner for OrbitSim test suite."""

import sys
import pytest

if __name__ == '__main__':
    args = sys.argv[1:] if len(sys.argv) > 1 else ['-v', 'tests']
    ret = pytest.main(args)
    sys.exit(ret)
