:Begin:
:Function: InitializeSolver
:Pattern: InitializeSolver[prob_String, n_Integer, s_Real, rL_Real, uL_Real, eL_Real, rR_Real, uR_Real, eR_Real]
:Arguments: {prob, n, s, rL, uL, eL, rR, uR, eR}
:ArgumentTypes: {String, Integer, Real, Real, Real, Real, Real, Real, Real}
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
