import aizynthfinder.aizynthfinder
from aizynthfinder.aizynthfinder import AiZynthFinder
import tensorflow as tf
import os

gpu_devices = tf.config.experimental.list_physical_devices('GPU')
for device in gpu_devices:
    tf.config.experimental.set_memory_growth(device, True)

if not os.path.exists("azf"):
    os.mkdir("azf")

# TODO: run the download public data script if azf file does not contain hdf5 files.

config = "azf/config.yml"
finder = AiZynthFinder(configfile=config)

finder.stock.select("zinc")
finder.expansion_policy.select("uspto")
finder.filter_policy.select("uspto")

finder.target_smiles = "Cc1cccc(c1N(CC(=O)Nc2ccc(cc2)c3ncon3)C(=O)C4CCS(=O)(=O)CC4)C"
finder.tree_search()

finder.build_routes()
stats = finder.extract_statistics()

print(stats)
# stats.get("feasibility")
