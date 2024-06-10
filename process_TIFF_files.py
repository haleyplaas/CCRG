# %%

from PIL import Image
import numpy as np 

tiff_data = Image.open('Q:\My Drive\Code Repositories\R\CCRG\SeaDAS files\CyanoIndices\L2022152.L3m_DAY_CYAN_CI_cyano_CYAN_CONUS_300m_8_3.tif')
Image.close()
cyan_index = np.array(tiff_data)

print(cyan_index)
# %%
