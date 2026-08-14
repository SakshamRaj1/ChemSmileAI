import h5py

with h5py.File("molecules_library.h5", "r") as f:

    def show_structure(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(
                f"{name:30s} "
                f"dtype={str(obj.dtype):10s} "
                f"shape={obj.shape}"
            )

    f.visititems(show_structure)

# import h5py

# h5_file = "molecules_library.h5"

# with h5py.File(h5_file, "r") as f:
#     print("Datasets:")
    
#     def inspect(name, obj):
#         if isinstance(obj, h5py.Dataset):
#             print(f"{name}")
#             print(f"  shape: {obj.shape}")
#             print(f"  dtype: {obj.dtype}")

#     f.visititems(inspect)