cd ../../build/Release

chmod +x gridgen

./gridgen ../../data/Figure22/grid_1.json ../../data/Figure22/config.json -t 0.008 -o "CSG" --tree ../../data/Figure22/doghead_800_int_tree.json
