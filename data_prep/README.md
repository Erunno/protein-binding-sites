# data_prep

This folder primarily provides abstractions of the Yu dataset and PDB files, and operations with them.

## Script [datasets_db.py](./datasets_db.py)

This script includes the following Python classes:
- [`ChainRecord`](./datasets_db.py): Represents one protein chain and offers many methods to shield the user from both the Yu dataset and computations like SASA and protrusion.
- [`LigandDataset`](./datasets_db.py): Gathers chains into datasets as defined by the Yu dataset.
- [`SeqDatasetDb`](./datasets_db.py): A class collecting all the datasets and providing methods to retrieve any `LigandDataset` you want.

## Script [pdb_files_db.py](./pdb_files_db.py)

This script focuses on the abstraction of 3D structures:
- [`Chain3dStructure`](./pdb_files_db.py): Provides a set of functions for protrusion, SASA, and closest neighbors, etc. This abstraction is primarily used by the [`ChainRecord`](./datasets_db.py).
- [`PdbFilesDb`](./datasets_db.py): Collects `Chain3dStructure`s and provides a set of functions for retrieving `Chain3dStructure`.

## File Cache

Due to the time-consuming nature of operations implemented by [`Chain3dStructure`](./pdb_files_db.py), we use a caching mechanism to avoid recomputation:
- [`FileCache`](./file_cache.py) is utilized with a `use_cache` function, allowing any object to be cached. See the example below for usage details.

```python
class BaseClass:
    def __init__(self):
        # class have to ready for caching
        # it has to define these two properties:
        #    __cache_key__: defines a unique id of the object
        #    __cashable__:  defines which methods are cashable
        self.__cache_key__ = 'key_to_cache'
        self.__cashable__ = [
            self.expensive_operation.__name__
        ]

    def expensive_operation(self, *args, **kwargs):
        print(f"Performing expensive operation... for {args}, {kwargs}")
        return f"Result of expensive operation ... for {args}, {kwargs}"
    

base_class = BaseClass()

# `use_cache` creates a proxy object `cached_base`
# that catches the call of all functions and search cache
# for the call with same parameters
with use_cache(base_class) as cached_base:

    # runs expensive operation
    res1 = cached_base.expensive_operation(1)
    print()
    
    # runs expensive operation
    res2 = cached_base.expensive_operation(1, test='test')
    print()
    
    # finds the result in cache and do not have ta call the 
    # expensive method
    res3 = cached_base.expensive_operation(1, test='test')
    print()

    print ('res 1:', res1)
    print ('res 2:', res2)
    print ('res 3:', res3)
```

We also defined the script  [`pdb_files_refresh_cache.py`](./pdb_files_refresh_cache.py) to prefill the cache.

## How to Use Data from Yu and PDB

Apart from the implementation of the classes discussed above we have also prepared multiple handy functions that takes the information you want (like embedings, prtrusion, sasa etc.) and transform it into the vectors that are directly usebale in machine learning so that no further preprocessing on the part of neural network implementation is needed. The typical usecase can be seen in the script [`network.v2.py`](../netws/network.v2.py) where it is used or see this example:

The folder also contains functions that transform information into vectors directly usable in machine learning, as demonstrated in the script [`network.v2.py`](../netws/network.v2.py) or the following example:


```python

# define the database
pdb_db = pdb_data.PdbFilesDb()
db = dataset.SeqDatasetDb()
db.set_pdb_db(pdb_db)

# load dataset for desired ligand
ds = db.get_dataset_for(ligand)

# obtain data you need
X_train, y_train, X_test, y_test = ds.get_train_test_data(
    [
        # these are utility functions that retrieves needed information
        # here our feature vectors gonna contain embeddings
        # from the embedder ESM, than protrusion for radius 8.5 A and lastly
        # the SASA value
        dataset.DataAccessors.embeddings("ESM"),
        dataset.DataAccessors.protrusion(8.5),
        dataset.DataAccessors.SASA_vector(),
    ],

    # also we need to filter chains that has valid 3D file
    # i.e. the PDB file defines the same amino acids sequence as the Yu dataset record
    filters=[dataset.Helpers.filter_chains_with_valid_3D_file]
)

# `X_train`, `X_test` are matrices containing the desired feature vectors 
#                     and each feature vector comprises embeddings, 
#                     protrusion for 8.5 and SASA
# `y_train`, `y_test` are label vectors determining ligandability
```

## Tests

We have created a series of tests for both the [`ChainRecord`](./datasets_db.py) and [`Chain3dStructure`](./pdb_files_db.py), found in the scripts [`pdb_files_tests.py`](./pdb_files_tests.py) and [`datasets_tests.py`](./datasets_tests.py).