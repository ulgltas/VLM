#!/usr/bin/env python3
# -*- coding: utf8 -*-
# Copyright 2019 Université de Liège
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


def create_workspace(name, origin):
    import os
    from distutils import dir_util

    cur_dir = os.getcwd()
    work_dir = name.replace(os.sep, "_")
    work_dir = work_dir.replace(".", "_")
    workspace = os.path.join(cur_dir, "workspace", work_dir)
    origin_models = os.path.join(cur_dir, origin, "models")
    models = os.path.join(workspace, "models")
    if not os.path.exists(models):
        os.makedirs(models)
    if os.path.exists(origin_models):
        dir_util.copy_tree(origin_models, models)
    os.chdir(workspace)
    return workspace

def add_path(name, usecurdir = False):
    import os, sys
    if usecurdir:   # Uses current path
        curdir = os.getcwd()
        path = os.path.join(curdir, name)
    else:           # Uses relative path from VLM directory
        vlm_dir = os.path.abspath(os.path.split(__file__)[0])
        path = os.path.join(vlm_dir, "..", "..", name)
    if os.path.isdir(path):
         print(("INFO: adding {} to PYTHONPATH".format(path)))
         sys.path.append(path)
    else:
        print(("INFO: {} not found!".format(path)))
