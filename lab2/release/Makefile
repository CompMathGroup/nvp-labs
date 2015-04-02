.PHONY: pack

pack:
	make -C ../build
	mkdir -p Windows_x86 Windows_x86-64 Unix_x86 Unix_x86-64
	cp ../gui/LabParabolic.nb Parabolic.nb
	cp ../build/win32/HeatLink.exe Windows_x86/
	cp ../build/win64/HeatLink.exe Windows_x86-64/
	cp ../build/unix32/HeatLink.exe Unix_x86/
	cp ../build/unix64/HeatLink.exe Unix_x86-64/
	zip -r Lab2.zip Parabolic.nb Unix_* Windows_*
