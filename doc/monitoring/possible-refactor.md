Data aggregation output and intervention deployment
===============

Data is potentially per-timestep, per-person, per mosquito category, per
treatment type, per other thing (grouping people into categories?). Some of
this data is output specific.

A very generic method: at every possible point a reportable piece of data is
available, call a user-defined function. This may temporarily aggregate data or
write it to an output array immediately. At the end of the time-step, another
user-defined function is called, which can, for example, be used to calculate
aggregates across people for the timestep (or other types of aggregates), using
arithmetic means or any other method. Is this fast enough, or would too much
data be generated? Would be interesting to test.

Output arrays: I suggest these are set up statically or in some cases with
dynamic rows (but statically defined dimensions). Output would be in the form
of multiple tables. Tables may include a single measure or multiple in every
cell. 2D tables may be output in a table-like format (current continuous
reporting), while 3D tables must be more generic (classic format, or maybe
something better compressed, if it reduces size even after zipping). Would then
have multiple output files. The current two files would each be a "table".

Cohorts and such: this system shouldn't be fully independent of intervention
deployment. Simplest thing would be to have multiple "cohorts" (or categories)
of people and assign each individual to a category at birth and potentially
change later, using user-defined functions. For higher complexity, could use
multiple labels instead of a single one, but this may not be necessary.

Intervention deployment is similar to cohort selection in that generally people
"have an intervention or not", except that the system remembers *when* the
intervention was administered. Could use user-defined functions to choose
intervention deployment at set times and ages based on available properties
(age, weight, recent infections, cohort, other interventions).
