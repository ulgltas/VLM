#! /usr/bin/env python
# -*- coding: utf8 -*-

def create_workspace(name, origin):
    import os
    from distutils import dir_util

    vlm_dir = os.path.abspath(os.path.split(__file__)[0])
    workspace = os.path.join(vlm_dir, "..", "workspace", name)
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