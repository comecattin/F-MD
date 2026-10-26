PREFIX ?= /usr/local
PYTHON ?= python3
COMMANDS = md_simulation fmd-plot-energies fmd-plot-trajectories

.PHONY: all install test clean
all:
	$(MAKE) -C src

install: all
	install -d "$(DESTDIR)$(PREFIX)/bin" "$(DESTDIR)$(PREFIX)/libexec/f-md"
	install -m 755 src/md_simulation "$(DESTDIR)$(PREFIX)/libexec/f-md/md_simulation"
	install -m 644 src/input_parser.py src/plot_energies.py src/plot_trajectories.py "$(DESTDIR)$(PREFIX)/libexec/f-md/"
	set -e; for command in $(COMMANDS); do \
		install -m 755 scripts/fmd-command "$(DESTDIR)$(PREFIX)/bin/$$command"; \
	done

test:
	bash tests/run.sh
	$(PYTHON) -m unittest discover -s tests -p 'test_*.py' -v

clean:
	$(MAKE) -C src clean
