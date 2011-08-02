Some of the documentation here is written in markup. While it should be very
readable as-is, it can also be converted to html.

Requirements: cmake, make, sed, pandoc

Method:
	mkdir build
	cd build
	ccmake .. (press c,c,g)
	make

The html files are generated in the same folder as the markdown files, while
the cmake build-stuff is out of the way in the build folder.
