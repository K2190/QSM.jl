CC = gcc  
CFLAGS = -std=c99 -fPIC -Wall -Wextra -O2 -g -DGAUSS_SEIDEL_RED_BLACK
LDFLAGS = -shared   
RM = rm -f   

TARGET_LIB = libjulia_qsm.a  

LIBS = -lm

SRCS = boundary_mask_l1_mex.c boundary_mask_mex.c coarsen_grid_mex.c \
       correct_mex.c fmg_mex.c gauss_seidel_mex.c gradb_adj_mex.c    \
       gradb_mex.c gradbm_adj_mex.c gradbm_mex.c gradc_mex.c gradcm_mex.c \
       gradf_adj_mex.c gradf_mex.c gradfm_adj_mex.c gradfm_mex.c \
       gradfp_adj_mex.c gradfp_mex.c lap1_mex.c lapmg_mex.c lapw_mex.c \
       mg_mex.c mgpcg_mex.c norm2_mex.c pd_tgv_mex.c pd_tv_mex.c prolong_mex.c \
       residual_mex.c restrict_mex.c mx_util.c mx_blas.c
       
OBJS = $(SRCS:.c=.o)

.PHONY: all
all: ${TARGET_LIB}

$(TARGET_LIB): $(OBJS) 
	$(CC) ${LDFLAGS} $(LIBS) -o $@ $^

$(SRCS:.c=.d):%.d:%.c
	$(CC) $(CFLAGS) -MM $< >$@

include $(SRCS:.c=.d)

.PHONY: clean
clean:
	-${RM} ${TARGET_LIB} ${OBJS} $(SRCS:.c=.d)
