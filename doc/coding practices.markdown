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


Data storage
------------------------

### Compile-time constants

Use `const static ty name = ...;` (at namespace level, `static` can be skipped).
See: https://stackoverflow.com/a/178259/314345

### Scenario constants

I.e. data which must be read from input but is constant throughout the
simulation. This data should be `static` and set by an `init` function.
It may need to use `std::unique_ptr<T>`.

The data "variable" should be private to its module.
If the data is only used within a single `.cpp` file, it should not be declared
in the header. If it is read from outside that file and is not performance
sensitive, it should also not be declared in the header and only accessed
indirectly through a function declared in the header. If reads are performance
sensitive, it should be declared in the header either via `extern` or
(preferably) as a private `static` struct/class member; either way reads should
be via an inline function defined in the header.

This type of data does not need checkpointing.

### Mutable simulation data

Most data falls into this category, and should be:

-   allocated on the stack
-   allocated within a struct/class which itself respects this policy
-   explicitly allocated on the heap by a smart pointer (e.g. `std::unique_ptr`)
    which is allocated according to this policy

The `Simulation` class is the root simulation object.
All data in this category likely needs checkpointing.

### Special mutable data

We make exceptions in a few cases:

-   simulation times, which are accessed *very* frequently and only updated
    at well defined times
-   the random number generator, which currently shares a single state between
    all users
-   reporting data, which is accessed from many different places

The last two will need attention for parallelisation of OM.


Forward declarations
--------------------------------

Forward declarations of classes, etc. just to avoid #including another file should only be used if:

* necessary to avoid a circular dependency of headers
* the header being written is included in a lot of compilation units (compilation performance)


Memory management
-----------------------------

There are two forms of memory management in C++:

*   do-it-yourself C-style with pointers
*   managed containers: vector, unique_ptr, etc.

I've tended towards the latter. Usually it's possible to arrange data into a
hierarchical structure, and in this case these containers will do most of the
work.

TODO: boost::shared_ptr has been used in the past to get around the limitations
of `std::auto_ptr`; now that we are using a later C++ standard we can migrate
away from some of these BOOST facilities.
