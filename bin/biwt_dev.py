import os
def biwt_dev_mode(biwt):
    BIWT_SPATIAL_DEV = os.getenv('BIWT_SPATIAL_DEV', 'False')
    BIWT_SPATIAL_DEV = BIWT_SPATIAL_DEV=='True'
    if BIWT_SPATIAL_DEV:
        self.column_line_edit.setText("cluster")
        # self.import_file("./data/visium_adata.h5ad")
        self.import_file("./data/Zhuang-ABCA-1-1.064_raw_wClusterAnnots.h5ad")
        # self.continue_from_import()
        # self.continue_from_spatial_query()
        # self.continue_from_edit()
        # self.window.process_window() # process rename window
        # self.set_cell_positions()
    else:
        self.import_file("./data/pbmc3k_clustered.h5ad")
        self.continue_from_edit()
        self.window.process_window() # process rename window
        self.window.process_window() # process cell count window
        # self.continue_from_rename()
        # self.set_cell_positions()