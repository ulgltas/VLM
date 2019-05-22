#!/usr/bin/env python
# -*- coding: utf8 -*-
# VLM run script

def create_workspace(name, origin):
    import os
    from distutils import dir_util

    vlm_dir = os.path.abspath(os.path.split(__file__)[0])
    workspace = os.path.join(vlm_dir, "workspace", name)
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
    path = os.path.join(vlm_dir, "..", name)
    if os.path.isdir(path):
         print("INFO: adding {} to PYTHONPATH".format(path))
         sys.path.append(path)
    else:
        print("INFO: {} not found!".format(path))

def main():
    add_path("geoGen")
    add_path("VLM/bin")

    import argparse, os

    parser = argparse.ArgumentParser()
    parser.add_argument("file", nargs=1, help="Python input file")
    parser.add_argument("--nogui", action="store_true", help="Disable GUI", default=False)
    args = parser.parse_args()


    vlm_dir = os.path.abspath(os.path.split(__file__)[0])

    names = args.file[0].split(os.sep)
    filename = names[-1]
    name = names[-2]
    origin = os.path.split(args.file[0])[0]

    add_path(os.path.join("VLM",origin))

    workspace = create_workspace(name, origin)
    print("Workspace: {}".format(workspace))

    

    filename = os.path.join(vlm_dir, "tests", name, filename)
    if os.path.isfile(filename):
        print("Executing file {}".format(filename))
        execfile(filename, globals(), locals())
    else:
        raise Exception("File not found: {}".format(filename))


if __name__ == "__main__":
    main()