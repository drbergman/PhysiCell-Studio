import os
def biwt_dev_mode(biwt):
    BIWT_SPATIAL_DEV = os.getenv('BIWT_SPATIAL_DEV', 'False')
    BIWT_SPATIAL_DEV = BIWT_SPATIAL_DEV=='True'
    if BIWT_SPATIAL_DEV:
        biwt.column_line_edit.setText("cluster")
        biwt.import_file("./data/visium_adata.h5ad")
        # biwt.continue_from_import()
        biwt.continue_from_spatial_query()
        biwt.continue_from_edit()
        biwt.window.process_window() # process rename window
        # biwt.set_cell_positions()
    else:
        biwt.import_file("./data/pbmc3k_clustered.h5ad")
        biwt.continue_from_edit()
        biwt.window.process_window() # process rename window
        biwt.window.process_window() # process cell count window
        # biwt.continue_from_rename()
        # biwt.set_cell_positions()