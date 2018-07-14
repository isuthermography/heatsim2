

all:
	python setup.py build

clean:
	python setup.py clean
	rm -rf build/
	rm -f *.pyc *~ */*~ */*.pyc */*/*~ */*/*.pyc doc/*~

commit: clean
	git add -A #hg addrem
	git commit -a #hg commit
