import pubchempy as pcp
from PIL import Image
import numpy as np

pcp.download('PNG', 'drug.png', 10096043, 'cid', overwrite=True)
img = Image.open('drug.png').convert('L')
img.save('./greyscale.png')
pixels = np.array(img)
pixels[pixels == 245] = 255
pixels[pixels < 245] = 0
pixels = pixels.flatten()

