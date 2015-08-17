// Barrier
//
//////////////////////////////////////////////////////////////////////


#ifndef _BARRIER_
#define _BARRIER_


#define PAD (4)
typedef struct{
volatile __declspec(align(4096)) unsigned _SIMPLE_BARRIER_TURN_[256*PAD];
volatile __declspec(align(4096)) unsigned _SIMPLE_BARRIER_FLAG_0;
volatile __declspec(align(4096)) unsigned _SIMPLE_BARRIER_COUNT_0[256];
volatile __declspec(align(4096)) unsigned _SIMPLE_BARRIER_FLAG_1;
volatile __declspec(align(4096)) unsigned _SIMPLE_BARRIER_COUNT_1[256];
} SIMPLE_BARRIER_TYPE;

inline void SIMPLE_BARRIER_INIT(SIMPLE_BARRIER_TYPE *SIMPLE_BARRIER)
{
    SIMPLE_BARRIER->_SIMPLE_BARRIER_FLAG_0 = 0;
    SIMPLE_BARRIER->_SIMPLE_BARRIER_FLAG_1 = 0;
    int i;
    for (i = 0; i < 256; i++)
    {
        SIMPLE_BARRIER->_SIMPLE_BARRIER_TURN_[i*PAD] = 0;
        SIMPLE_BARRIER->_SIMPLE_BARRIER_COUNT_0[i] = 0;
        SIMPLE_BARRIER->_SIMPLE_BARRIER_COUNT_1[i] = 0;
    }
}

inline void SIMPLE_BARRIER(SIMPLE_BARRIER_TYPE *SIMPLE_BARRIER, int ID, int MAX)
{
		int i;
    if (SIMPLE_BARRIER->_SIMPLE_BARRIER_TURN_[ID*PAD] == 0)
    {
        SIMPLE_BARRIER->_SIMPLE_BARRIER_TURN_[ID*PAD] = 1;
        if (ID == 0)
        {
            for (i = 1; i < MAX; i++)
            {
                while(SIMPLE_BARRIER->_SIMPLE_BARRIER_COUNT_0[i] == 0);
                SIMPLE_BARRIER->_SIMPLE_BARRIER_COUNT_1[i] = 0;
            }
            SIMPLE_BARRIER->_SIMPLE_BARRIER_FLAG_1 = 0;
            SIMPLE_BARRIER->_SIMPLE_BARRIER_FLAG_0 = 1;
        }
        else
        {
            SIMPLE_BARRIER->_SIMPLE_BARRIER_COUNT_0[ID] = 1;
            while(SIMPLE_BARRIER->_SIMPLE_BARRIER_FLAG_0 == 0);
        }
    }
    else
    {
        SIMPLE_BARRIER->_SIMPLE_BARRIER_TURN_[ID*PAD] = 0;
        if (ID == 0)
        {
            for (i = 1; i < MAX; i++)
            {
                while(SIMPLE_BARRIER->_SIMPLE_BARRIER_COUNT_1[i] == 0);
                SIMPLE_BARRIER->_SIMPLE_BARRIER_COUNT_0[i] = 0;
            }
            SIMPLE_BARRIER->_SIMPLE_BARRIER_FLAG_0 = 0;
            SIMPLE_BARRIER->_SIMPLE_BARRIER_FLAG_1 = 1;
        }
        else
        {
            SIMPLE_BARRIER->_SIMPLE_BARRIER_COUNT_1[ID] = 1;
            while(SIMPLE_BARRIER->_SIMPLE_BARRIER_FLAG_1 == 0);
        }
    }
}

#define PAD (4)
#define GROUP (4)
typedef struct {
volatile __declspec(align(4096)) unsigned _TREE_BARRIER_TURN_[256*PAD];
volatile __declspec(align(4096)) unsigned _TREE_BARRIER_FLAG_0[4][256*PAD];
volatile __declspec(align(4096)) unsigned _TREE_BARRIER_COUNT_0[4][256*PAD];
volatile __declspec(align(4096)) unsigned _TREE_BARRIER_FLAG_1[4][256*PAD];
volatile __declspec(align(4096)) unsigned _TREE_BARRIER_COUNT_1[4][256*PAD];
} TREE_BARRIER_TYPE;

inline void TREE_BARRIER_INIT(TREE_BARRIER_TYPE *TREE_BARRIER)
{
		int i;
    for (i = 0; i < 256; i++)
    {
        TREE_BARRIER->_TREE_BARRIER_TURN_[i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_FLAG_0[0][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_FLAG_0[1][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_FLAG_0[2][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_FLAG_0[3][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_COUNT_0[0][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_COUNT_0[1][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_COUNT_0[2][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_COUNT_0[3][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_FLAG_1[0][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_FLAG_1[1][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_FLAG_1[2][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_FLAG_1[3][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_COUNT_1[0][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_COUNT_1[1][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_COUNT_1[2][i*PAD] = 0;
        TREE_BARRIER->_TREE_BARRIER_COUNT_1[3][i*PAD] = 0;
    }
}

inline void TREE_BARRIER(TREE_BARRIER_TYPE *TREE_BARRIER, int ID, int MAX)
{
		int SET;
    if (TREE_BARRIER->_TREE_BARRIER_TURN_[ID*PAD] == 0)
    {
        TREE_BARRIER->_TREE_BARRIER_TURN_[ID*PAD] = 1;
        int LEVEL = 0;
        int STRIDE = 1;
        int MASK = GROUP - 1;
        while (1)
        {
            int WAITER = ID & MASK;
            int SUBSET = ID & (~MASK);
            if (WAITER)
            {
                TREE_BARRIER->_TREE_BARRIER_COUNT_0[LEVEL][ID*PAD] = 1;
                while(TREE_BARRIER->_TREE_BARRIER_FLAG_0[LEVEL][SUBSET*PAD] == 0);
                while (LEVEL > 0)
                {
                    LEVEL--;
                    TREE_BARRIER->_TREE_BARRIER_FLAG_1[LEVEL][ID*PAD] = 0;
                    TREE_BARRIER->_TREE_BARRIER_FLAG_0[LEVEL][ID*PAD] = 1;
                }
                break;
            }
            else
            {
                int CHILD = ID + STRIDE;
                for (SET = 1; (SET < GROUP) && (CHILD < MAX); SET++)
                {
                    while(TREE_BARRIER->_TREE_BARRIER_COUNT_0[LEVEL][CHILD*PAD] == 0);
                    TREE_BARRIER->_TREE_BARRIER_COUNT_1[LEVEL][CHILD*PAD] = 0;
                    CHILD += STRIDE;
                }
            }

            LEVEL++;
            STRIDE *= GROUP;
            MASK *= GROUP;
            if (STRIDE >= MAX)
            {
                while (LEVEL > 0)
                {
                    LEVEL--;
                    TREE_BARRIER->_TREE_BARRIER_FLAG_1[LEVEL][ID*PAD] = 0;
                    TREE_BARRIER->_TREE_BARRIER_FLAG_0[LEVEL][ID*PAD] = 1;
                }
                break;
            }
        }
    }
    else
    {
        TREE_BARRIER->_TREE_BARRIER_TURN_[ID*PAD] = 0;
        int LEVEL = 0;
        int STRIDE = 1;
        int MASK = GROUP - 1;
        while (1)
        {
            int WAITER = ID & MASK;
            int SUBSET = ID & (~MASK);
            if (WAITER)
            {
                TREE_BARRIER->_TREE_BARRIER_COUNT_1[LEVEL][ID*PAD] = 1;
                while(TREE_BARRIER->_TREE_BARRIER_FLAG_1[LEVEL][SUBSET*PAD] == 0);
                while (LEVEL > 0)
                {
                    LEVEL--;
                    TREE_BARRIER->_TREE_BARRIER_FLAG_0[LEVEL][ID*PAD] = 0;
                    TREE_BARRIER->_TREE_BARRIER_FLAG_1[LEVEL][ID*PAD] = 1;
                }
                break;
            }
            else
            {
                int CHILD = ID + STRIDE;
                for (SET = 1; (SET < GROUP) && (CHILD < MAX); SET++)
                {
                    while(TREE_BARRIER->_TREE_BARRIER_COUNT_1[LEVEL][CHILD*PAD] == 0);
                    TREE_BARRIER->_TREE_BARRIER_COUNT_0[LEVEL][CHILD*PAD] = 0;
                    CHILD += STRIDE;
                }
            }

            LEVEL++;
            STRIDE *= GROUP;
            MASK *= GROUP;
            if (STRIDE >= MAX)
            {
                while (LEVEL > 0)
                {
                    LEVEL--;
                    TREE_BARRIER->_TREE_BARRIER_FLAG_0[LEVEL][ID*PAD] = 0;
                    TREE_BARRIER->_TREE_BARRIER_FLAG_1[LEVEL][ID*PAD] = 1;
                }
                break;
            }
        }
    }
}


#endif
