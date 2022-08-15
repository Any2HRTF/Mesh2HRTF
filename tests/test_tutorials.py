import os
import utils

# define and check paths to your Blender versions
blender_path = utils.blender_paths(2)[0]
ValueError(blender_path)

base_dir = os.path.dirname(__file__)

# utils.instal_blender_addons()