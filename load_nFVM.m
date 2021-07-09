global RUN_OPTIMISATIONS
global VECTORIZE

Sqrt = @(field) fields.sqrt(field);
Det = @(field) fields.det(field);
Log = @(field) fields.log(field);
Inv = @(field) fields.inv(field);
Ones = @(field) fields.ones(field);

RUN_OPTIMISATIONS = false;
VECTORIZE = false;

loadPaths
