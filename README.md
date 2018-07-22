# OrbitSim
## My Vision
This is currently a 2D orbit plotter on Python. 
I hope to make this into a space program simulator game, where the user "programs" each probe. As opposed to simulating the craft in real time, like KSP does, time will go in steps as specified by the user. Burns will be single events. Once the user has programmed the mission plan, they can incorporate it into the game. The program will go through simulate each active mission plan until a particular time that the user specifies. 

I'd like the style to be minimalist, where a player can easily follow what's happening. If a feature adds complexity, I don't want it. For example, rather than designing a spacecraft by parts, I'm thinking of treating the spacecraft as a single object with resources and maybe cost. This will lessen the learning curve and encourage player imagination.

Abstraction is also important. I want the future user to have access to simple tools with which they can design their mission. I want this to be simple when possible but allow the user to have the standard freedoms of programming. On the other hand, I don't want the user to be able to just type "craft.orbit.periapsis = 25.23456e3" and have it be so. But they might want to know that information for their own program, maybe calculating how much to burn if they want apoapsis to be some value. Or another user might want to burn until apoapsis is some value without learning how to calculate that.

## Contributing
I'd like to try opening this up for other people to contribute to. It's one of the few time's I'm interested in the final product and not just the programming process. I'm hoping it'll be a game worth sharing and for real people to try out.
I'm new at using GitHub, so bear with me while I figure out how to collaborate.

## Log 1: July 2018
Currently this is a simple 2D Orbit Simulator. You can specify an orbit, calculate orbit elements, and most importantly calculate the state at a particular time. Alternately, you can specify a state at a particular time, and it'll calculate the orbit that it'll travel. 

It was quite an ordeal getting the math to work for weird cases. Many sources tell you how to do it for the simple elliptical orbit case, but I also wanted flybys and orbiting counter-clockwise, and I wanted it done simply. That is why I wanted it in 2D in the first place.

I imagine the next step will be defining a basic craft that can change orbits. We need to allow a craft in a hyperbolic trajectory away from Earth to then orbit the Sun. Then it should check if it's eventually influenced by some other planet. There's probably some simple way to do this. Of course we're using circles of influence. 

For now, a craft can be an object with an orbit. A burn can be a change in velocity. Fuel can be infinite, but it should be programmed with that feature in mind.
