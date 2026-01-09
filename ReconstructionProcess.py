import napari
import tifffile as tiff
from tkinter import Tk, filedialog

# 1. 打开文件选择对话框
root = Tk()
root.withdraw()  # 不显示主窗口

file_path = filedialog.askopenfilename(
    title="Select a 3D TIFF / OME-TIFF file",
    filetypes=[("TIFF files", "*.tif *.tiff *.ome.tif")]
)

if not file_path:
    raise SystemExit("No file selected.")

# 2. 读取数据
volume = tiff.imread(file_path)
print("Data shape:", volume.shape)

# 3. Napari 3D 显示
viewer = napari.Viewer(ndisplay=3)
viewer.add_image(
    volume,
    name="Zebrafish",
    rendering="mip"  # 可改为 'iso' / 'translucent'
)
napari.run()
