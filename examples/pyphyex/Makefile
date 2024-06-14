clean:
	rm -rf firstSedimTest CB02_GAUS __pycache__ test_microphy_cold firstShallowTest firstTurbTest

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

all: ice_adjust microphy sedimentation shallow turb
