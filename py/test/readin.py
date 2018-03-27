num_atoms = len(open('coords.dat', 'r').readlines())
coords = np.zeros((num_atoms, 2))

fn = open('coords.dat', 'r')
for i in range(num_atoms):
    c_info = fn.readline().split()
    print c_info
fn.close()