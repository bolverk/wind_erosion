def main():

    import matplotlib
    matplotlib.use('Agg')
    import pylab
    import glob
    import re
    import numpy
    import h5py
    import tqdm
    import imageio

    def extract_snapshot_number(fpath):

        return int(re.search(r'_(\d+).',fpath).group(1))

    ss_files = sorted(glob.glob('output/snapshot_*.h5'),
                      key=extract_snapshot_number)
    for fname in tqdm.tqdm(ss_files):
        with h5py.File(fname,'r') as f:
            pylab.tricontourf(f['geometry']['x_coordinate'],
                              f['geometry']['y_coordinate'],
                              numpy.log10(f['hydrodynamic']['density']),
                              50)
            pylab.title('t = %10.2f' % numpy.array(f['time'])[0])
            pylab.colorbar()
            pylab.axis('equal')
            pylab.savefig(fname.replace('.h5','.png'))
            pylab.close()

    with imageio.get_writer('density.gif',
                            mode='I',
                            duration=0.5) as writer:
        for fname in tqdm.tqdm(ss_files):
            image = imageio.imread(fname.replace('.h5','.png'))
            writer.append_data(image)

if __name__ == '__main__':

    main()
