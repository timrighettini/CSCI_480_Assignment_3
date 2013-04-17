======================
==========
Tim Righttini's Assignment #3 Readme: CSCI_480 Computer Graphics (CG)
==========
======================

------
Hello:
------

This Program contains a raytracer that can render spheres and triangles, specifically.

The program completes the following tasks, EXCLUDING extra credit:
1.  Uniformly send out rays from the camera location.  Final images are 640x480, and the fov is 60 degrees.
2.  Write the intersection code.  This has been done for triangles and spheres.
3.  Implement the illumination equations for triangle and sphere shading
4.  Create still images showing off your ray tracer.  These still are located within this folder.

Extra Credit Includes the following:
1.  Implemented antialiasing that is adjustable within the code.  If you want to change the amount of extra samples taken per pixel, just modify the "sampleNumber" value located at line 89 of assign3.cpp.  If the sample number is two, 2 x 2 sub pixels will be averaged in the end result for anti-aliasing.
2.  Implemented soft shadowing that is also adjustable within the code.  If you want to change the amount of random lights located around any given light source, change the parameter "numRandomLights" located at line 90 in assign3.cpp.
3.  Initially, sampleNumber is set to 1 (and cannot be set below 1, else no rays will draw to the screen) and numRandomLights is set to 0.  To see the effects of 1. and 2. outside the stills, you will have to modify the parameters aforementioned in their proper places.

---------
Controls:
---------

None: Just watch the camera roll!

-------------
Special Notes
-------------

*Notes about stills:
Images 000.jpg through 004.jpg only show the base functionality of the raytracer.
Image 005.jpg shows off the soft shadows functionality of my raytracer with test2, with 100 random lights located at all of the light sources.
Image 006.jpg shows off the soft shadows and anti-aliasing functionality of my raytracer with test2, with 100 random lights and sampleNumber set to 4.  Notice how crisp the image looks.
Image 007.jpg used the same parameters are image 006.jpg, but with the spheres scene.
Image 008.jpg used 10 lights per light source and the sampleNumber set to 2, but with the table scene.  This still has a stronger focus on soft shadows versus anti-aliasing.
Image 009.jpg shows off the soft shadows and anti-aliasing functionality of my raytracer with test2, with 25 random lights and sampleNumber set to 5.  Notice how crisp the image looks.

*Note: Trying to use a sampleNumber of 8 or more will lead to a compile error due to memory concerns.  You may be able to use 7 and get away with it, but 36 pixels samples per pixel is more than enough.  I also would have liked to do a special render with the letters, but it was just taking too long to render.  I have the base render as  for the letters image 000.jpg, though, to show that I could render it regularly without any EC, but it is possible to render it with EC parameters.

*Note: The screenfile is initially set to test1.scene, but I have txt files of the other scenes within my submission.  TO test these, just copy and paste the content to the screenfile.

Besides for the notes and what I have mentioned above, this is everything you need to know about the program.

Thanks! 

-------
The End
-------
