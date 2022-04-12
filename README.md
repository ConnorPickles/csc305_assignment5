# Assignment 5

## Environment

- **OS:** Ubuntu 20.04.3 LTS (running in WSL 2)
- **Compiler:** Clang 10.0.0

## Implementation

I implemented all required parts of this assignment. There are too many results to practically put in the readme, so please find them in the `img_results` folder. There are separate sub-folders for the orthogonal and perspective results.

I noticed the lighting in my results were slightly different than the lighting in the provided results. My best guess for the cause of this is the specular exponent; when I set it to something like 16 instead of 256, the noticeable light reflections disappeared (mainly from the pv shading) and my results looked more like the provided ones. Either way, the lighting is very similar and you can see in my results that the lighting changes as the bunny rotates.