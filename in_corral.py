#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

if __name__ == "__main__":
    os.environ.setdefault("CORRAL_SETTINGS_MODULE", "rbogus.settings")

    from corral import cli
    cli.run_from_command_line()
