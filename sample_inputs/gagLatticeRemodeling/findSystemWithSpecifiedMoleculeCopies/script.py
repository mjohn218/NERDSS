import numpy as np
import os
import shutil


copies = 2500

all_trajs = []
all_trajs_pol = []

path = "./restart67p"

isExist = os.path.exists(path)

if not isExist:
    os.makedirs(path)

for i in range(1, 241):
    one_traj_time = np.genfromtxt('./%d/sim_time.dat' % (i))
    one_traj = np.genfromtxt('./%d/sim_gagMem.dat' % (i))
    one_traj_pol = np.genfromtxt('./%d/sim_polMem.dat' % (i))
    for oneindex in range(0, len(one_traj_time)):
        if one_traj[oneindex] >= copies:
            checkpoint = round(one_traj_time[oneindex], 2)
            for twoindex in range(0, len(one_traj_time)):
                if abs(one_traj_time[twoindex] - checkpoint) < 1E-20:
                    #print(i, checkpoint, one_traj[twoindex], one_traj_pol[twoindex], one_traj[twoindex] +
                    #      one_traj_pol[twoindex], one_traj_pol[twoindex]/(one_traj[twoindex]))
                    j = int(checkpoint * 1000) * 2000
                    path = f"./restart67p/{i}"
                    isExist = os.path.exists(path)
                    if not isExist:
                        os.makedirs(path)

                    src = f"./{i}/restart{j}.dat"
                    dst = f"./restart67p/{i}/restart.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"../{i}/rng_state"
                    dst = f"./restart67p/{i}/rng_state"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_time.dat"
                    dst = f"./restart67p/{i}/sim_time.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_gagMem.dat"
                    dst = f"./restart67p/{i}/sim_gagMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_maxMem.dat"
                    dst = f"./restart67p/{i}/sim_maxMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_polMem.dat"
                    dst = f"./restart67p/{i}/sim_polMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_fragment.dat"
                    dst = f"./restart67p/{i}/sim_fragment.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
            break

copies = 1250

all_trajs = []
all_trajs_pol = []

path = "./restart33p"

isExist = os.path.exists(path)

if not isExist:
    os.makedirs(path)

for i in range(1, 241):
    one_traj_time = np.genfromtxt('./%d/sim_time.dat' % (i))
    one_traj = np.genfromtxt('./%d/sim_gagMem.dat' % (i))
    one_traj_pol = np.genfromtxt('./%d/sim_polMem.dat' % (i))
    for oneindex in range(0, len(one_traj_time)):
        if one_traj[oneindex] >= copies:
            checkpoint = round(one_traj_time[oneindex], 2)
            for twoindex in range(0, len(one_traj_time)):
                if abs(one_traj_time[twoindex] - checkpoint) < 1E-20:
                    #print(i, checkpoint, one_traj[twoindex], one_traj_pol[twoindex], one_traj[twoindex] +
                    #      one_traj_pol[twoindex], one_traj_pol[twoindex]/(one_traj[twoindex]))
                    j = int(checkpoint * 1000) * 2000
                    path = f"./restart33p/{i}"
                    isExist = os.path.exists(path)
                    if not isExist:
                        os.makedirs(path)

                    src = f"./{i}/restart{j}.dat"
                    dst = f"./restart33p/{i}/restart.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"../{i}/rng_state"
                    dst = f"./restart33p/{i}/rng_state"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_time.dat"
                    dst = f"./restart33p/{i}/sim_time.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_gagMem.dat"
                    dst = f"./restart33p/{i}/sim_gagMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_maxMem.dat"
                    dst = f"./restart33p/{i}/sim_maxMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_polMem.dat"
                    dst = f"./restart33p/{i}/sim_polMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_fragment.dat"
                    dst = f"./restart33p/{i}/sim_fragment.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
            break
        
copies = 625

all_trajs = []
all_trajs_pol = []

path = "./restart17p"

isExist = os.path.exists(path)

if not isExist:
    os.makedirs(path)

for i in range(1, 241):
    one_traj_time = np.genfromtxt('./%d/sim_time.dat' % (i))
    one_traj = np.genfromtxt('./%d/sim_gagMem.dat' % (i))
    one_traj_pol = np.genfromtxt('./%d/sim_polMem.dat' % (i))
    for oneindex in range(0, len(one_traj_time)):
        if one_traj[oneindex] >= copies:
            checkpoint = round(one_traj_time[oneindex], 2)
            for twoindex in range(0, len(one_traj_time)):
                if abs(one_traj_time[twoindex] - checkpoint) < 1E-20:
                    #print(i, checkpoint, one_traj[twoindex], one_traj_pol[twoindex], one_traj[twoindex] +
                    #      one_traj_pol[twoindex], one_traj_pol[twoindex]/(one_traj[twoindex]))
                    j = int(checkpoint * 1000) * 2000
                    path = f"./restart17p/{i}"
                    isExist = os.path.exists(path)
                    if not isExist:
                        os.makedirs(path)

                    src = f"./{i}/restart{j}.dat"
                    dst = f"./restart17p/{i}/restart.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"../{i}/rng_state"
                    dst = f"./restart17p/{i}/rng_state"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_time.dat"
                    dst = f"./restart17p/{i}/sim_time.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_gagMem.dat"
                    dst = f"./restart17p/{i}/sim_gagMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_maxMem.dat"
                    dst = f"./restart17p/{i}/sim_maxMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_polMem.dat"
                    dst = f"./restart17p/{i}/sim_polMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_fragment.dat"
                    dst = f"./restart17p/{i}/sim_fragment.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
            break
        
copies = 1875

all_trajs = []
all_trajs_pol = []

path = "./restart50p"

isExist = os.path.exists(path)

if not isExist:
    os.makedirs(path)

for i in range(1, 241):
    one_traj_time = np.genfromtxt('./%d/sim_time.dat' % (i))
    one_traj = np.genfromtxt('./%d/sim_gagMem.dat' % (i))
    one_traj_pol = np.genfromtxt('./%d/sim_polMem.dat' % (i))
    for oneindex in range(0, len(one_traj_time)):
        if one_traj[oneindex] >= copies:
            checkpoint = round(one_traj_time[oneindex], 2)
            for twoindex in range(0, len(one_traj_time)):
                if abs(one_traj_time[twoindex] - checkpoint) < 1E-20:
                    #print(i, checkpoint, one_traj[twoindex], one_traj_pol[twoindex], one_traj[twoindex] +
                    #      one_traj_pol[twoindex], one_traj_pol[twoindex]/(one_traj[twoindex]))
                    j = int(checkpoint * 1000) * 2000
                    path = f"./restart50p/{i}"
                    isExist = os.path.exists(path)
                    if not isExist:
                        os.makedirs(path)

                    src = f"./{i}/restart{j}.dat"
                    dst = f"./restart50p/{i}/restart.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"../{i}/rng_state"
                    dst = f"./restart50p/{i}/rng_state"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_time.dat"
                    dst = f"./restart50p/{i}/sim_time.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_gagMem.dat"
                    dst = f"./restart50p/{i}/sim_gagMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_maxMem.dat"
                    dst = f"./restart50p/{i}/sim_maxMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_polMem.dat"
                    dst = f"./restart50p/{i}/sim_polMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_fragment.dat"
                    dst = f"./restart50p/{i}/sim_fragment.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
            break
        
copies = 3125

all_trajs = []
all_trajs_pol = []

path = "./restart84p"

isExist = os.path.exists(path)

if not isExist:
    os.makedirs(path)

for i in range(1, 241):
    one_traj_time = np.genfromtxt('./%d/sim_time.dat' % (i))
    one_traj = np.genfromtxt('./%d/sim_gagMem.dat' % (i))
    one_traj_pol = np.genfromtxt('./%d/sim_polMem.dat' % (i))
    for oneindex in range(0, len(one_traj_time)):
        if one_traj[oneindex] >= copies:
            checkpoint = round(one_traj_time[oneindex], 2)
            for twoindex in range(0, len(one_traj_time)):
                if abs(one_traj_time[twoindex] - checkpoint) < 1E-20:
                    #print(i, checkpoint, one_traj[twoindex], one_traj_pol[twoindex], one_traj[twoindex] +
                    #      one_traj_pol[twoindex], one_traj_pol[twoindex]/(one_traj[twoindex]))
                    j = int(checkpoint * 1000) * 2000
                    path = f"./restart84p/{i}"
                    isExist = os.path.exists(path)
                    if not isExist:
                        os.makedirs(path)

                    src = f"./{i}/restart{j}.dat"
                    dst = f"./restart84p/{i}/restart.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"../{i}/rng_state"
                    dst = f"./restart84p/{i}/rng_state"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_time.dat"
                    dst = f"./restart84p/{i}/sim_time.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_gagMem.dat"
                    dst = f"./restart84p/{i}/sim_gagMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_maxMem.dat"
                    dst = f"./restart84p/{i}/sim_maxMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_polMem.dat"
                    dst = f"./restart84p/{i}/sim_polMem.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
                    src = f"./{i}/sim_fragment.dat"
                    dst = f"./restart84p/{i}/sim_fragment.dat"
                    try:
                        shutil.copyfile(src, dst)
                    except:
                        print(f"error in copyfile: {src} to {dst}")
            break