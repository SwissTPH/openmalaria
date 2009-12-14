Coding practices
===========


File naming
------------------

If a file contains a single class, the name excluding extension should correspond exactly to the
class name (currently abbreviations sometimes used; avoid in future).

Use the file extensions .cpp (for compilation units) and .hpp (for include files) (most include
files currently use .h however).


Namespaces
-------------------

Aim to implement the following namespace structure:

    namespace OM { namespace SUBDIR { namespace MODULE {
	...
    } } }

OM is the namespace for the OpenMalaria model. SUBDIR and MODULE should correspond exactly to the
sub-directory of include/model (assuming one sub-dir deep) and the file name excluding extension.

(This is new and largely unimplemented.)
