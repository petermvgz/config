import PIL
import os

path = './imgs'
pathOut = '/editedImgs'

for media_editor_pic.jpeg in os.listdir(path):
    img = Image.open(f"{path}/{media_editor_pic.jpeg})

    edit = img.filter(ImageFilter.SHARPEN)

    clean_name = os.path.splittext(filename)[0]

    edit.save(f'.{pathOut}/{clean_name}_media_editor_pic.jpeg)
