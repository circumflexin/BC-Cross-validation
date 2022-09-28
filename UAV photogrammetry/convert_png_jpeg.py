import os

from PIL import Image


n = 1
for file in os.listdir(".\\"):
    if file.endswith(".png"):
        im = Image.open(file)
        rgb_im = im.convert('RGB')
        rgb_im.save("%s.jpg" % file[:-4])
        n += 1