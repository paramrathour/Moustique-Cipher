#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('display', 'latex')


# In[2]:


from sage.matrix.berlekamp_massey import berlekamp_massey
import numpy as np
import bisect
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


# In[3]:


def init(key, F, register_lengths, row_count, r_1_num):
    CCSR = list(random_matrix(F, 1, register_lengths[0]+1).row(0))    # 1-indexed
    r_1 = random_matrix(F, r_1_num, register_lengths[1]+3)
    r_2 = list(random_matrix(F, 1, register_lengths[2]).row(0))
    r_3 = list(random_matrix(F, 1, register_lengths[3]).row(0))
    
    sum = np.zeros(len(row_count), dtype = int)
    sum[0] = row_count[0]
    for i in range(1, len(row_count)):
        sum[i] = row_count[i] + sum[i-1]
    
    mapping = {}
    mapping[0,0] = 0
    column_count = np.zeros(len(key)+1, dtype = int)
    column_count[0] = 1
    for i in range(1,register_lengths[0]+1):
        idx = bisect.bisect_left(sum, i)
        j = 48 + i - sum[idx]
        mapping[idx, j] = i
        column_count[j] += 1
    for i in range(2 + column_count[len(key)]):
        j = len(key) - i
        if (1, j) not in mapping:
            mapping[1, j] = 0
    
    return CCSR, r_1, r_2, r_3, mapping, column_count


# In[4]:


def g(i, a, b, c, d):
    if i == 0:
        return a + b + c + d
    elif i == 1:
        return 1 + a + b + c * (d + 1)
    elif i == 2:
        return a * (1 + b) + c * (1 + d)
    else:
        return 0


# In[5]:


def parameters(i,j):
    value = (j-i) % 6
    if value == 1 or value == 4:
        return 0, 2*(j-i-1)/3, j-2
    elif value == 2 or value == 5:
        return 1, j-4, j-2
    elif value == 3:
        return 1, 0, j-2
    else:
        return 1, j-5, 0


# In[6]:


def state_update(F, key, CCSR, r_1, r_2, r_3, mapping, column_count, register_lengths, r_1_num, row_count):
    CCSR_new = list(random_matrix(F, 1, register_lengths[0]+1).row(0))    # 1-indexed
    for j in range(1, len(key)+1):
        for i in range(column_count[j]):
            x, v, w = parameters(i, j)
            if j < 3:
                v = 0
                w = 0
                c = 0
                d = 0
            a = CCSR[mapping[i%column_count[j-1], j-1]]
            b = key[j-1]
            c = CCSR[mapping[i%column_count[v], v]]
            d = CCSR[mapping[i%column_count[w], w]]
            if j == len(key):
                x = 2
                b = CCSR[mapping[0, j-1-i]]
                c = CCSR[mapping[i%column_count[j-2], j-2]]
                d = CCSR[mapping[1, j-2-i]]
            CCSR_new[mapping[i, j]] = g(x, a, b, c, d)
    CCSR = CCSR_new;
    for i in range(register_lengths[1]):
        a = CCSR[register_lengths[0] - i]
        b = CCSR[i + int(register_lengths[0] * 18/128)]
        c = CCSR[int(register_lengths[0] * 113/128) - i]
        d = CCSR[i + 1]
        r_1[0, (4*i) % register_lengths[1]] = g(1, a, b, c, d)
    for j in range(1, r_1_num):
        for i in range(register_lengths[1]):
            a = r_1[j-1, i]
            b = r_1[j-1, i+3]
            c = r_1[j-1, i+1]
            d = r_1[j-1, i+2]
            r_1[j, (4*i) % register_lengths[1]] = g(1, a, b, c, d)
    for i in range(register_lengths[2]):
        a = r_1[r_1_num-1, 4*i]
        b = r_1[r_1_num-1, 4*i+3]
        c = r_1[r_1_num-1, 4*i+1]
        d = r_1[r_1_num-1, 4*i+2]
        r_2[i] = g(1, a, b, c, d)
    for i in range(register_lengths[3]):
        a = r_2[4*i]
        b = r_2[4*i+1]
        c = r_2[4*i+2]
        d = r_2[4*i+3]
        r_3[i] = g(0, a, b, c, d)
    return CCSR, r_1, r_2, r_3


# In[7]:


def output_stream(key, CCSR, r_1, r_2, r_3):
    return sum(r_3)


# In[8]:


def generate_stream(F, key, IV, CCSR, r_1, r_2, r_3, mapping, column_count, register_lengths, r_1_num, row_count, size):
    key_combined = [key[i] + IV[i] for i in range(len(key))]
    output = []
    for i in range(min(size, len(key_combined))):
        output.append(key_combined[i])
    for i in range(len(key_combined), size):
        CCSR, r_1, r_2, r_3 = state_update(F, key_combined, CCSR, r_1, r_2, r_3, mapping, column_count, register_lengths, r_1_num, row_count)
        output.append(output_stream(key, CCSR, r_1, r_2, r_3))    
    return output


# In[9]:


def generate_data(IV_num, size):
    F = GF(2)
    register_lengths = [64, 53, 12, 3]
    r_1_num = 2
    row_count = [48, 6 , 3, 3, 1, 1, 1, 1]
    output = []
    IVs = []
    key_length = 48
    key = random_matrix(F, 1, key_length)
    key = list(key.row(0))
    for i in range(IV_num):
        CCSR, r_1, r_2, r_3, mapping, column_count = init(key, F, register_lengths, row_count, r_1_num)
        IV = list(random_matrix(F, 1, key_length).row(0))
        IVs.append(IV)
        output.append(generate_stream(F, key, IV, CCSR, r_1, r_2, r_3, mapping, column_count, register_lengths, r_1_num, row_count, size))
    return key, IVs, output


# In[10]:


set_random_seed(0)
key, IVs, output = generate_data(100, 1024)


# In[11]:


linear_complexities = []
for o in output:
    linear_complexities.append(berlekamp_massey(o+o).degree())


# In[12]:


fig, ax = plt.subplots()
ax.plot(np.arange(1, len(linear_complexities)+1), linear_complexities, '-om')
ax.set_xlabel('Run Number', fontsize = 20)
ax.set_ylabel('Linear Complexity', fontsize = 20)
fig.set_size_inches(15, 10, forward=True)
ax.set_ylim([1018, 1025])
ax.set_title('Variation of Linear Complexity with Runs', fontsize = 25)
plt.savefig('1.pdf')


# In[13]:


bins_num = len(np.unique(linear_complexities))
fig, ax = plt.subplots()
ax.hist(linear_complexities, bins = bins_num, rwidth = 0.8, color = 'violet')
ax.set_xlabel('Linear Complexity', fontsize = 20)
ax.set_ylabel('Count', fontsize = 20)
fig.set_size_inches(15, 10, forward=True)
ax.set_xlim([1018, 1025])
ax.set_title('Histogram of Linear Complexity Values', fontsize = 25)
plt.savefig('2.pdf')


# In[14]:


r = ZZ.random_element(0,99)
linear_complexities_one = []
for i in range(2, len(output[r])):
    linear_complexities_one.append(berlekamp_massey(o[1:i]+o[1:i]).degree())


# In[15]:


fig, ax = plt.subplots()
ax.plot(linear_complexities_one, '-m')
ax.set_xlabel('Stream Length', fontsize = 20)
ax.set_ylabel('Linear Complexity', fontsize = 20)
fig.set_size_inches(15, 10, forward=True)
ax.set_title('Linear Complexity as Stream Length Increases for a Random Run {}'.format(r), fontsize = 25)
plt.savefig('3.pdf')

