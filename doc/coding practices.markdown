Coding practices
===========


Naming conventions
-------------------------------

* constants — USE_CAPS
* variables & functions — useLowercase
* types — UseCamelCase


File naming
------------------

If a file contains a single class, the name excluding extension should correspond exactly to the
class name (currently abbreviations sometimes used; avoid in future).

Use the file extensions .cpp (for compilation units) and .hpp (for include files) (most include
files currently use .h however).


Namespaces
-------------------

Aim to implement the following namespace structure:

    namespace OM { namespace SUBDIR {
	...
    } }

or, when it's preferable to have a separate namespace from the parent dir,

    namespace OM { namespace SUBDIR { namespace MODULE {
	...
    } } }


OM is the namespace for the OpenMalaria model. SUBDIR and MODULE should correspond exactly to the
sub-directory of include/model (assuming one sub-dir deep) and the file name excluding extension.

(This is new and largely unimplemented.)


Static members
------------------------

Avoid making any data requiring checkpointing static, so it's clear all static data doesn't change
after program initialization. One exception to this is a little of the data in Global, which is
really only there to avoid passing the values constantly.


Forward declarations
--------------------------------

Forward declarations of classes, etc. just to avoid #including another file should only be used if:

* necessary to avoid a circular dependency of headers
* the header being written is included in a lot of compilation units (compilation performance)
