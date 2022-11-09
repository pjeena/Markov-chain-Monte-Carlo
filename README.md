# Markov-chain-Monte-Carlo

Everyone knows about a magnet 🧲 right. We know that a magnet has a north pole and a south pole. A north(south) pole attracts a south(north) pole and repels a north(south) pole. Let's suppose a magnet to be an arrow ⬆️ where the top of arrow (∧) is the north pole and bottom is south pole. We know that
the same poles repel each other and opposite poles attract(**thus the magnets like to stay  ⬆️ ⬇️ with each other** , let's call this **fundamental happiness condition**). Going one step ahead, we but one more magnet
from the market and we want to place them on a triangle(such that they live happily after). How would one do that ??

**TRIAL** :: We place the first  magnet(⬆️) on one corner of the triangle, the second magnet(⬇️) on the second corner, what should be the orientation of the third magnet? Is it \textbf{$\uparrow$} or \textbf{$\downarrow$}. 

![frus](https://user-images.githubusercontent.com/56345075/200874408-d3e9443f-cf41-4d52-815b-435dd2b790ab.png)


Case (i) :: If it is ⬆️ ?  or Case (ii) :: If it is ⬇️ ?  Both the possibilities are explained in the figure.


Case (i)            |  Case (ii)
:-------------------------:|:-------------------------:
![case11](https://user-images.githubusercontent.com/56345075/200870882-2f8c63b6-d7c3-4676-842b-a45412e67084.png)  |  ![case22](https://user-images.githubusercontent.com/56345075/200870937-3ddd5622-dea9-4c6d-9a5f-7cd74ad8952d.png)

But, we see that in both the cases one pair of magnets(1 and 3 in case (i); 2 and 3 in case (ii)) is not satisying the **fundamental happiness condition**; 
⬆️ ⬇️. This is called FRUSTRATION in physics where the magnets cannot live happily with each other for a given geometric shape; obviously the name comes 
from the english vocabulary. Moreover this situation can arise at any of the three corners of a triangle thus giving 2 * 3 =  6 mathematical combinations shown below :

![Ising-spins-on-triangular-lattice-Six-possible-configurations-for-ground-state-are](https://user-images.githubusercontent.com/56345075/200875431-6e66e739-aefc-425b-a163-6e109e2810da.png)

Sometimes, a little bit of frustration makes life interesting 😂. Although controversial, this statement is certainly correct in physics, where frustration has opened an entire field of research with materials of exotic properties with potential applications to develop essential technologies for quantum computers. One of these technologies is “topologically protected” qubits. Qubits are information processors, but these qubits are easily perturbed by outside influences. To make a topologically protected qubit, we could use the collective motion of our frustrated electrons to encode information in a way that does not care much about the environment. One can think of it like how a knot remains tied regardless of what might be tickling the ends of the string.

In my thesis work, we go beyond a triangle and put electons(tiny magnets with an intrinsic spin) on a tetrahedron based geometry shown below. 
![id18089](https://user-images.githubusercontent.com/56345075/200860318-6918c70c-7ac4-4003-8823-d0b15775b998.jpg)

The number of combinations/configrations possible for this type of system is huge and intractable analytically. Our aim is to calculate the average values of random variables(like energy etc) for our system thus we need to know the probability/likelihood of each configration. The expected value is calculated by multiplying each of the possible outcomes by the likelihood each outcome will occur and then summing all of those values. This is insanely difficult in our case since we can't track these many outcomes. Therefore, we need to use a very popular clever sampling technique called as MARKOV CHAIN MONTE CARLO.







