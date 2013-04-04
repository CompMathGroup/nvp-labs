:Begin:
:Function: InitializeSolver
:Pattern: InitializeSolver[solver_String, k_Real, c_Real, prob_String, n_Integer, s_Real]
:Arguments: {solver, k, c, prob, n, s}
:ArgumentTypes: {String, Real, Real, String, Integer, Real}
:ReturnType: Real
:End:

:Begin:
:Function: AddScheme
:Pattern: AddScheme[deltas_List]
:Arguments: {deltas}
:ArgumentTypes: {Real64List}
:ReturnType: Integer
:End:

:Begin:
:Function: DoSteps
:Pattern: DoSteps[count_Integer]
:Arguments: {count}
:ArgumentTypes: {Integer}
:ReturnType: Manual
:End:

int main(int argc, char **argv) {
	return MLMain(argc, argv);
}
