INCLUDES = -I,+oUnit/_build,-I,/Users/farr/Documents/code/ocaml_randomize/_build,-I,+gsl,-I,+diff/_build

CFLAGS = -cflags $(INCLUDES),-g
LFLAGS = -lflags $(INCLUDES),-g

LIBS = unix,diff,oUnit,randomize,bigarray,gsl

OCB_DIRS = -Is analysis

OCB = ocamlbuild $(CFLAGS) $(LFLAGS) $(OCB_DIRS)
OCB_EXE = $(OCB) -libs $(LIBS) $(OCB_DIRS)

.PHONY: all
all: lib

.PHONY: lib
lib:
	$(OCB) niv.cma niv.cmxa

.PHONY: sho-dEs
sho-dEs:
	$(OCB_EXE) sho_dEs.native

.PHONY: kep-dEs
kep-dEs:
	$(OCB_EXE) kep_dEs.native

.PHONY: ecc-kep-dEs
ecc-kep-dEs:
	$(OCB_EXE) ecc_kep_dEs.native

.PHONY: bh-dEs
bh-dEs:
	$(OCB_EXE) bh_dEs.native

.PHONY: bh-traj-err
bh-traj-err:
	$(OCB_EXE) bh_traj_err.native

.PHONY: clean
clean:
	$(OCB) -clean

.PHONY: test
test:
	$(OCB_EXE) run_tests.native --