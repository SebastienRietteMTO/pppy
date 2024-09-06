clean:
	rm -rf firstSedimTest CB02_GAUS __pycache__ test_microphy_cold firstShallowTest firstTurbTest \
		LIMA test_LIMA_cold

ice_adjust:
	./comp_ice_adjust.py

microphy:
	./comp_microphy.py

sedimentation:
	./comp_sedimentation.py

shallow:
	./comp_shallow.py

turb:
	./comp_turb.py

lima_adjust:
	./comp_lima_adjust.py

lima:
	./comp_lima.py

all: ice_adjust microphy sedimentation shallow turb lima_adjust lima
