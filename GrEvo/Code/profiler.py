import pstats
#python3 -m cProfile -o profile.txt GrEvolution.py ../Data/G_Data/Brindabellaspis.txt ../Data/G_Data/Meemannia.txt ../Data/C_Data/Brindabellaspis.txt ../Data/C_Data/Meemannia.txt 0.1
with open("Profile.txt", "w") as f:
    ps = pstats.Stats("SYPA_Profile.txt", stream=f)
    ps.sort_stats('cumulative')
    ps.print_stats()