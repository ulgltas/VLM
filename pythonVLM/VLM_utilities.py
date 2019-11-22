#! /usr/bin/env python
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

    temp_dir = os.path.join(os.path.split(__file__)[0], "..")
    vlm_dir = os.path.abspath(temp_dir)
    work_dir = origin.replace(os.sep, "_")
    workspace = os.path.join(vlm_dir, "workspace", work_dir)
    origin_models = os.path.join(vlm_dir, origin, "models")
    models = os.path.join(workspace, "models")
    if not os.path.exists(models):
        os.makedirs(models)
    if os.path.exists(origin_models):
        dir_util.copy_tree(origin_models, models)
    os.chdir(workspace)
    return workspace

def add_path(name):
    import os, sys
    vlm_dir = os.path.abspath(os.path.split(__file__)[0])
    path = os.path.join(vlm_dir, "..", "..", name)
    if os.path.isdir(path):
         print("INFO: adding {} to PYTHONPATH".format(path))
         sys.path.append(path)
    else:
        print("INFO: {} not found!".format(path))
