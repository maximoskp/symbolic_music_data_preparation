#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 17:43:46 2019

@author: maximoskaliakatsos-papakostas
"""

import score2np as s2n
import os
cwd = os.getcwd()

folderName = cwd + '/bc404'
parts_for_surface = 'all'
time_res = 16
ts_num = 4
ts_den = 4

m = s2n.get_time_sig_parts_np_from_folder(folderName, parts_for_surface, time_res, ts_num, ts_den, print_progress=True, transpose=True, range_trim=True)