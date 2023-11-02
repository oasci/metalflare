# 03-data

TODO:

## Organization

The organization is loose, but specific templates for your data types may exist.

Each experiment/dataset gets one subfolder.
Here comes data that should be saved for long-term storage, published, and/or shared with collaborators.
The experiment subfolder starts with a three-digit number.

Derivative data that can be easily recreated (and should not be stored for the long term) should go in a different folder, 990_processed data.
I suggest splitting the data between experiments there, too.

## Contents

```{toctree}
:maxdepth: 1
:glob:

*
*/README
```
