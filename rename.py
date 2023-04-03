#%%

folder_path = r'C:\Users\jiaxiny\Documents\stanford_work\del'
old_ext = ".txt.md"
new_ext = ".txt"

import os

def convert_extension(folder_path, old_ext, new_ext):
    for file_name in os.listdir(folder_path):
        if file_name.endswith(old_ext):
            old_file_path = os.path.join(folder_path, file_name)
            new_file_path = os.path.join(folder_path, file_name[:-len(old_ext)] + new_ext)
            os.rename(old_file_path, new_file_path)
            print(f"Renamed {old_file_path} to {new_file_path}")

# Example usage
convert_extension(folder_path, old_ext, new_ext)

# %%
