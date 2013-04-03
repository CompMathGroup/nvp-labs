.PHONY: pack

pack:
	make -C ../build
	mkdir -p Windows_x86 Windows_x86-64 Unix_x86 Unix_x86-64
	cp ../gui/LabHyperbolic.nb Hyperbolic.nb
	cp ../build/win32/GasLink.exe Windows_x86/
	cp ../build/win64/GasLink.exe Windows_x86-64/
	cp ../build/unix32/GasLink.exe Unix_x86/
	cp ../build/unix64/GasLink.exe Unix_x86-64/
	zip -r Lab1.zip Hyperbolic.nb Unix_* Windows_*
