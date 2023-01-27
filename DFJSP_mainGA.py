# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 01:02:50 2022

@author: parkh
"""

import random
import copy
import pandas as pd
import numpy as np
from datetime import datetime
from collections import Counter
import plotly.express as px
from datetime import datetime
import matplotlib.pyplot as plt
params = {
    'MUT': 1,  # 변이확률(%)
    'END' : 0.9,  # 설정한 비율만큼 chromosome이 수렴하면 탐색을 멈추게 하는 파라미터 (%)
    'POP_SIZE' : 100,  # population size 10 ~ 100
    'RANGE' : 10, # chromosome의 표현 범위, 만약 10이라면 00000 00000 ~ 11111 11111까지임
    'NUM_OFFSPRING' : 5, # 한 세대에 발생하는 자식 chromosome의 수
    'SELECTION_PRESSURE' : 3, # 선택연산의 선택압
    'job_seq' : [y for x in range(1,5) for y in range(1,11)],
    'factory_seq' : [1,2]
    # 원하는 파라미터는 여기에 삽입할 것
    }
class JAYA_FJSP():
    def __init__(self, parameters):
        self.eps = 0.5
        self.eps2 = 0.1
        self.eps3 = 0.2
        self.params = {}
        for key, value in parameters.items():
            self.params[key] = value
        self.p_table = pd.read_csv('C:/Users/parkh/DFJSP3.csv',index_col=(0)) #job과 operation을 기록한 테이블
        self.s_table = pd.read_csv('C:/Users/parkh/DFJSP_SETUP.csv', index_col=(0)) #setup time 테이블
        self.d_table = pd.read_csv('C:/Users/parkh/DFJSP_del2.csv', index_col=(0))
        self.job_endTime={'j1':0, 'j2':0, 'j3':0, 'j4':0, 'j5':0, 'j6':0,'j7':0, 'j8':0, 'j9':0, 'j10':0} # job의 끝나는 지점을 등록
        self.machine_endTime={'M1':0,'M2':0,'M3':0,'M4':0, 'M5':0,'M6':0} # machine의 끝나는 지점을 등록
        self.machine_prejob={'M1':"j0", 'M2':"j0",'M3':"j0", 'M4':"j0",'M5':"j0", 'M6':"j0",'M7':"j0", 'M8':"j0"}
        self.job_preOperation={'1':1,'2':1,'3':1,'4':1,'5':1,'6':1,'7':1,'8':1,'9':1,'10':1}
        self.job_of_factory={'1':['M1','M2','M3'], '2':['M4','M5','M6']}
        self.job_max_op = {'1':4, '2':4, '3':4, '4':4, '5':4, '6':4, '7':4, '8':4, '9':4, '0':4}
        
    def reset(self):
        self.job_endTime={'j1':0, 'j2':0, 'j3':0, 'j4':0, 'j5':0, 'j6':0,'j7':0, 'j8':0, 'j9':0, 'j10':0} # job의 끝나는 지점을 등록
        self.machine_endTime={'M1':0,'M2':0,'M3':0,'M4':0, 'M5':0,'M6':0,'M7':0,'M8':0} # machine의 끝나는 지점을 등록
        self.machine_prejob={'M1':"j0", 'M2':"j0",'M3':"j0", 'M4':"j0",'M5':"j0", 'M6':"j0",'M7':"j0", 'M8':"j0"}
        self.job_preOperation={'1':1,'2':1,'3':1,'4':1,'5':1,'6':1,'7':1,'8':1,'9':1, '10':1}
    def gannt_chart(self, population):
        """
        assignment = [[1,'M3'],[1,'M2'],[10,'M1'],[9,'M3'],[9,'M3'],[7,'M1'],[10,'M2'],
                      [3,'M3'],[9,'M2'],[10,'M1'],[7,'M3'],[1,'M3'],[3,'M2'],[7,'M1'],
                      [10,'M2'],[9,'M3'],[3,'M1'],[3,'M2'],[1,'M1'],[7,'M3'],[2,'M5'],
                      [6,'M6'],[5,'M4'],[8,'M6'],[2,'M6'],[6,'M5'],[8,'M5'],[6,'M4'],
                      [6,'M4'],[4,'M6'],[4,'M6'],[5,'M5'],[2,'M4'],[8,'M6'],[4,'M5'],
                      [2,'M6'],[5,'M4'],[5,'M4'],[8,'M5'],[4,'M6']]
        
        assignment = [[3,'M3'],[5,'M1'],[3,'M2'],[6,'M3'],[5,'M2'],[7,'M1'],[4,'M3'],
                      [4,'M3'],[6,'M2'],[3,'M1'],[7,'M3'],[4,'M2'],[5,'M1'],[6,'M3'],
                      [7,'M1'],[5,'M2'],[4,'M3'],[7,'M3'],[3,'M2'],[6,'M1'],[8,'M6'],
                      [10,'M4'],[1,'M6'],[10,'M5'],[2,'M4'],[9,'M6'],[9,'M6'],[1,'M5'],
                      [10,'M4'],[2,'M6'],[8,'M5'],[8,'M4'],[9,'M5'],[1,'M6'],[1,'M5'],
                      [2,'M4'],[9,'M6'],[2,'M6'],[8,'M5'],[10,'M4']]
        
        assignment = [[2,"M2"],[7,"M5"],[4,"M3"],[6,"M6"],[1,"M1"],[5,"M1"],[8,"M3"],
                      [3,"M6"],[10,"M4"],[1,"M2"],[10,"M5"],[5,"M2"],[6,"M5"],[2,"M3"],
                      [1,"M3"],[7,"M6"],[4,"M3"],[3,"M5"],[7,"M4"],[8,"M2"],[10,"M4"],
                      [4,"M2"],[9,"M6"],[2,"M1"],[8,"M3"],[9,"M6"],[8,"M2"],[9,"M5"],
                      [3,"M6"],[3,"M5"],[6,"M4"],[2,"M3"],[5,"M1"],[7,"M6"],[4,"M3"],
                      [1,"M1"],[5,"M1"],[6,"M4"],[10,"M4"],[9,"M6"]]
        """
        assignment = []
        for j in range(40):
            machine = population[2][j]
            factory = self.machine_factory(machine[-1])
            assignment.append([population[0][j],population[2][j]])
        print(assignment)
        
        plotlydf = pd.DataFrame([],columns=['Task','Start','Finish','Resource']) #간트차트로 보여주기 위한 데이터프레임
        i=0 #간트차트의 인덱싱을 위한 숫자
        j=0
        for job_num,machine in assignment:  #['j11','M2']의 형태에서 잡과 머신을 가져옴
            job = 'j'+str(job_num)        #'j11'의 형태를 j1로 
            job_op = job+str(self.job_preOperation[str(job_num)])
            self.job_preOperation[str(job_num)] += 1
            df2_sorted = self.s_table[job] #셋업테이블에서 job에 해당하는 컬럼을 가져옴
            setup_time=df2_sorted.loc[self.machine_prejob[machine]] #컬럼에서 machine에 세팅되어있던 job에서 변경유무 확인
            time = max(self.machine_endTime[machine] ,self.job_endTime[job]) #machine과 job의 순서 제약조건을 지키기 위해 더 큰 값을 설정함
            df_sorted = self.p_table[machine] #p_time테이블에서 현재 machine에 해당하는 열을 가져옴
            p_time = df_sorted.loc[job_op] #해당하는 job과 operation의 시간을 가져옴
            start = datetime.fromtimestamp(time*3600) #포매팅 해줌
            time = time+setup_time # 프로세스타임과 셋업타임을 더해줌
            p_start=datetime.fromtimestamp(time*3600)
            time = time+p_time
            end = datetime.fromtimestamp(time*3600) #끝나는 시간 포매팅
            plotlydf.loc[j] = dict(Task=job, Start=p_start, Finish=end, Resource=machine) #간트차트를 위한 딕셔너리 생성, 데이터프레임에 집어넣음
            if setup_time !=0:
                j+=1
                plotlydf.loc[j] = dict(Task="setup", Start=start, Finish=p_start, Resource=machine) #간트차트를 위한 딕셔너리 생성, 데이터프레임에 집어넣음
            j += 1
            i += 1 #데이터 프레임 인덱싱 증가
            self.machine_endTime[machine]=time #기계의 끝나는 시간 설정
            self.job_endTime[job]=time #job의 끝나는 시간 설정
            self.machine_prejob[machine] = job #현재 어떤 machine에서 어떤 job을 수행했는지 기록      
        self.reset()
        plotlydf2 = plotlydf.sort_values(by=['Resource','Task'], ascending=False)
        fig = px.timeline(plotlydf2, x_start="Start", x_end="Finish", y="Resource", color="Task", width=1000, height=400)
        fig.show()
    def get_fittness(self,scheduling_seq,routing_seq):
        time_list=[]
        for j in range(40):
            job_number = scheduling_seq[j]
            machine = routing_seq[j]
            
            job = "j"+str(job_number)
            jobOp=job+str(self.job_preOperation[str(job_number)])
            self.job_preOperation[str(job_number)] += 1
            
            setup_list = self.s_table[job]
            setup_time=setup_list.loc[self.machine_prejob[machine]]
            start_time = max(self.machine_endTime[machine] ,self.job_endTime[job])
            p_list = self.p_table[machine]
            p_time = p_list.loc[jobOp]
            end_time = start_time +setup_time+p_time
            time_list.append([start_time,end_time])
            self.machine_endTime[machine]=end_time #기계의 끝나는 시간 설정
            self.job_endTime[job]=end_time #job의 끝나는 시간 설정
            self.machine_prejob[machine] = job #현재 어떤 machine에서 어떤 job을 수행했는지 기록
        all_values = self.machine_endTime.values()
        c_max=max(all_values)
        all_values = list(all_values)
        k=0
        critical_machine=""
        for i in range(6):
            if all_values[i] > k:
                k=all_values[i]
                critical_machine = str(i+1)
        return c_max, critical_machine
    def anneal_eps(self):
        self.eps -=0.001
        self.eps2 +=0.001
        self.eps2 = min(self.eps2, 0.6)
        self.eps3 +=0.002
        self.eps3 = min(self.eps3, 0.8)
    def Algorithm2(self, time_list, c_max):
        critical_operation=[]
        k=len(time_list)
        for i in range(k):
            if time_list[k-1-i][2] == c_max:
                critical_operation.append(time_list[k-1-i][0])
                c_max = time_list[k-1-i][1]
        critical_operation.sort()
        return critical_operation    
    def init_OS(self):
        random.shuffle(self.params['job_seq'])
        scheduling_seq=copy.deepcopy(self.params['job_seq'])
        return scheduling_seq
    def init_MA_Random(self, scheduling_seq, f_routing_seq): # 이게 걍 랜덤하게 초기화
        routing_seq=[]
        for operation in scheduling_seq:
            operation2 = "j" + str(operation) + str(self.job_preOperation[str(operation)])
            self.job_preOperation[str(operation)] += 1
            random.shuffle(self.job_of_factory[str(f_routing_seq[operation-1])])
            for i in self.job_of_factory[str(f_routing_seq[operation-1])]:
                a = self.p_table[i].loc[operation2]
                if a != 0:
                    routing_seq.append(i)
                    break
        return routing_seq
    def init_MA_LS(self, scheduling_seq, f_routing_seq): #이건 SPT초기화
        routing_seq=[]
        for operation in scheduling_seq:
            operation2 = "j" + str(operation) + str(self.job_preOperation[str(operation)])
            self.job_preOperation[str(operation)] += 1
            machine = self.least_time_machine(operation,operation2,f_routing_seq)
            routing_seq.append(machine)
        self.reset()
        return routing_seq
    def Local_Machine_Routing(self, scheduling_seq, routing_seq ,f_routing_seq, job): #이게 local 초기화
        routing_seq2=[]
        job_list = []
        f = f_routing_seq[job-1]
        for i in range(len(f_routing_seq)):
            if f_routing_seq[i] == f:
                job_list.append(i+1)
        for i in range(len(scheduling_seq)):
            operation = scheduling_seq[i]
            operation2 = "j" + str(operation) + str(self.job_preOperation[str(operation)])
            self.job_preOperation[str(operation)] += 1
            if operation not in job_list:
                machine = routing_seq[i]
                time = max(self.machine_endTime[machine] ,self.job_endTime["j"+str(job)])
                df_sorted = self.p_table[machine]
                p_time = df_sorted.loc[operation2]
                df2_sorted = self.s_table[operation2[:2]]
                setup_time=df2_sorted.loc[self.machine_prejob[machine]]
                start = time
                end = start+setup_time+p_time
                self.machine_endTime[machine]=end #기계의 끝나는 시간 설정
                self.job_endTime[job]=end #job의 끝나는 시간 설정
                self.machine_prejob[machine] = "j"+str(job) #현재 어떤 machine에서 어떤 job을 수행했는지 기록
            else:
                machine = self.least_time_machine(operation,operation2,f_routing_seq)
            routing_seq2.append(machine)
        self.reset()
        return routing_seq2
    def Local_Machine_Routing_Swap(self, scheduling_seq, routing_seq ,f_routing_seq, job, job2): #이건 factory swap 떄문에
        routing_seq2=[]
        job_list = []
        f = f_routing_seq[job-1]
        f2 = f_routing_seq[job2-1]
        for i in range(len(f_routing_seq)):
            if f_routing_seq[i] == f2:
                job_list.append(i+1)
        for i in range(len(f_routing_seq)):
            if f_routing_seq[i] == f:
                job_list.append(i+1)
        for i in range(len(scheduling_seq)):
            operation = scheduling_seq[i]
            operation2 = "j" + str(operation) + str(self.job_preOperation[str(operation)])
            self.job_preOperation[str(operation)] += 1
            if operation not in job_list:
                machine = routing_seq[i]
                time = max(self.machine_endTime[machine] ,self.job_endTime["j"+str(job)])
                df_sorted = self.p_table[machine]
                p_time = df_sorted.loc[operation2]
                df2_sorted = self.s_table[operation2[:2]]
                setup_time=df2_sorted.loc[self.machine_prejob[machine]]
                start = time
                end = start+setup_time+p_time
                self.machine_endTime[machine]=end #기계의 끝나는 시간 설정
                self.job_endTime[job]=end #job의 끝나는 시간 설정
                self.machine_prejob[machine] = "j"+str(job) #현재 어떤 machine에서 어떤 job을 수행했는지 기록
            else:
                machine = self.least_time_machine(operation,operation2,f_routing_seq)
            routing_seq2.append(machine)
        self.reset()
        return routing_seq2
    def least_time_machine(self, job,operation2,f_routing_seq): #그 시점에서 가장 작은 machine 선택해주는거
        best_machine = ""
        max_endTime=10000
        machine_seq = self.job_of_factory[str(f_routing_seq[job-1])]
        for i in machine_seq:
            time = max(self.machine_endTime[i] ,self.job_endTime["j"+str(job)])
            df_sorted = self.p_table[i]
            p_time = df_sorted.loc[operation2]
            if p_time != 0:
                df2_sorted = self.s_table[operation2[:2]]
                setup_time=df2_sorted.loc[self.machine_prejob[i]]
                start = time
                end = start+setup_time+p_time
                if max_endTime>end:
                    max_endTime = end
                    best_machine = i
        self.machine_endTime[best_machine]=max_endTime #기계의 끝나는 시간 설정
        self.job_endTime[job]=max_endTime #job의 끝나는 시간 설정
        self.machine_prejob[best_machine] = "j"+str(job) #현재 어떤 machine에서 어떤 job을 수행했는지 기록
        return best_machine
    
    def Two_Point_OS_LS(self, solution): #LS 세가지
        coin = random.randint(0,2)
        scheduling_seq,factory_seq, routing_seq, fittness, critical_machine = solution
        s_v = copy.deepcopy(scheduling_seq)
        r_v = copy.deepcopy(routing_seq)
        operation_seq = [i for i in range(40)]
        k = random.choice(operation_seq)
        a = operation_seq.pop(k)
        b = random.choice(operation_seq)
        if a > b:
            a,b = b,a
        s_v2 = scheduling_seq[a:b+1]
        r_v2 = routing_seq[a:b+1]
        s_v3 = copy.deepcopy(s_v2)
        if coin == 0:
            s_v3.reverse()
            r_v3 = []
            for i in range(b-a+1):
                for j in range(b-a+1):
                    if s_v3[i] == s_v2[j]:
                        r_v3.append(r_v2[j])
                        s_v2[j] = 0
                        break
            k=0
            for i in range(a,b+1):
                s_v[i] = s_v3[k]
                r_v[i] = r_v3[k]
                k+=1
            c_max, critical_machine2 = self.get_fittness(s_v, r_v)
            solution = [s_v,factory_seq ,r_v, c_max, critical_machine2]
            return solution
        elif coin == 1:
            s_v3[0] , s_v3[-1] = s_v3[-1], s_v3[0]
            r_v3 = []
            for i in range(b-a+1):
                for j in range(b-a+1):
                    if s_v3[i] == s_v2[j]:
                        r_v3.append(r_v2[j])
                        s_v2[j] = 0
                        break
            k=0
            for i in range(a,b+1):
                s_v[i] = s_v3[k]
                r_v[i] = r_v3[k]
                k+=1
            c_max, critical_machine2 = self.get_fittness(s_v, r_v)
            solution = [s_v,factory_seq ,r_v, c_max, critical_machine2]
            return solution
        elif coin == 2:
            k = s_v3.pop(-1)
            s_v3.insert(0,k)
            r_v3 = []
            for i in range(b-a+1):
                for j in range(b-a+1):
                    if s_v3[i] == s_v2[j]:
                        r_v3.append(r_v2[j])
                        s_v2[j] = 0
                        break
            k=0
            for i in range(a,b+1):
                s_v[i] = s_v3[k]
                r_v[i] = r_v3[k]
                k+=1
            c_max, critical_machine2 = self.get_fittness(s_v, r_v)
            solution = [s_v,factory_seq ,r_v, c_max, critical_machine2]
            return solution
    def MOX_operator(self,dad_ch2, mom_ch2): #유일한 교차연산
        mom_ch = copy.deepcopy(mom_ch2)
        dad_ch = copy.deepcopy(dad_ch2)
        point1 = random.randint(0, 39)
        point2 = random.randint(point1, 40)
        dad_list = []
        offspring = [-1 for i in range(40)]
        offspring2 = [-1 for i in range(40)]
        for i in range(point1,point2): #리스트를 뽑아
            dad_list.append(dad_ch[0][i])
        for i in range(len(dad_list)): # 인덱싱을 찾아서 -1로 바꿔줘
            for j in range(len(mom_ch[0])):
                if mom_ch[0][j] == dad_list[i]:
                    mom_ch[0][j] = -1
                    break
        for i in dad_list: #그 인덱싱에다가 집어넣어
            for j in range(40):
                if mom_ch[0][j] == -1:
                    offspring[j] = i
                    mom_ch[0][j] = 0
                    break
        for i in range(40): # 나머지를 엄마에서 집어넣어
            if offspring[i] == -1:
                offspring[i] = mom_ch[0][i]
        for i in range(40): # 여기는 할당 따라가기
            for j in range(40):
                if offspring[i] == dad_ch[0][j]:
                    offspring2[i] = dad_ch[2][j]
                    dad_ch[0][j] = -1
                    break
        c_max, critical_machine = self.get_fittness(offspring, offspring2)
        self.reset()
        off_cho = [offspring, dad_ch[1], offspring2, c_max, critical_machine]
        return off_cho
        
    def init_FA(self): #공장 초기화
        f_routing_seq=[]
        for i in range(10):
            f_routing_seq.append(random.choice(self.params['factory_seq']))
        return f_routing_seq
    def Local_FA_Move_LS(self, offspring):
        scheduling_seq ,factory_seq ,routing_seq, fittness, critical_machine = offspring
        a = Counter(factory_seq)
        factory_seq2 = copy.deepcopy(factory_seq)
        min_list = [key for key,value in a.items() if min(a.values()) == value]
        machine = "M"+str(critical_machine)
        j_list =[]
        for i in range(40): 
            if routing_seq[i] == machine:
                j_list.append(scheduling_seq[i])
        job = random.choice(j_list)
        factory_seq2[job-1] = min_list[0]
        routing_seq2 = self.Local_Machine_Routing(scheduling_seq, routing_seq, factory_seq2, job)
        self.reset()
        fittness2, critical_machine2 = self.get_fittness(scheduling_seq, routing_seq2)
        solution = [scheduling_seq, factory_seq2, routing_seq2, fittness2, critical_machine2]
        self.reset()
        return solution
    def machine_factory(self, machine): #해당 머신이 어느 공장인지 반환
        factory = ''
        if machine == '1' or machine == '2' or machine == '3':
            factory = '1'
        elif machine == '4' or machine == '5' or machine == '6':
            factory = '2'
        else:
            factory = '3'
        return factory
    def Local_FA_Swap_LS(self, offspring): #바뀐 공장 두개만 local로 초기화
        scheduling_seq ,factory_seq ,routing_seq, fittness, critical_machine = offspring
        factory_seq2 = copy.deepcopy(factory_seq)
        machine = "M"+str(critical_machine)
        j_list =[]
        not_j_list = []
        for i in range(40): 
            if routing_seq[i] == machine:
                j_list.append(scheduling_seq[i])
            else:
                not_j_list.append(scheduling_seq[i])
        job = random.choice(j_list)
        job2 = random.choice(not_j_list)
        factory_seq2[job-1], factory_seq2[job2-1] = factory_seq[job2-1], factory_seq[job-1]
        routing_seq2 = self.Local_Machine_Routing_Swap(scheduling_seq, routing_seq, factory_seq2, job,job2)
        self.reset()
        fittness2, critical_machine2 = self.get_fittness(scheduling_seq, routing_seq2)
        solution = [scheduling_seq, factory_seq2, routing_seq2, fittness2, critical_machine2]
        self.reset()
        return solution
    def SPT_MA_LS(self, offspring): # 과부화 걸린 기계에서 p_time이 짧은거로 바꿈
        coin = random.randint(1, 1)
        scheduling_seq ,factory_seq ,routing_seq, fittness, critical_machine, critical_operation = offspring
        if coin == 1:
            critical_machine_index = []
            random.shuffle(critical_operation)
            factory = self.machine_factory(critical_machine)    
            stop = False
            for i in critical_machine_index:
                k=0
                for j in range(i+1):
                    if scheduling_seq[j] == scheduling_seq[i]:
                        k+=1
                job_op = 'j'+str(scheduling_seq[i])+str(k)
                machine = routing_seq[i]
                for j in self.job_of_factory[factory]:
                    if j != machine and self.p_table[j].loc[job_op] < self.p_table[machine].loc[job_op] and self.p_table[j].loc[job_op] != 0 :
                        routing_seq[i] = j
                        stop = True
                if stop:
                    break
            fittness, critical_machine,critical_operation = self.get_fittness(scheduling_seq, routing_seq)
            self.reset()
            solution = [scheduling_seq, factory_seq, routing_seq, fittness, critical_machine,critical_operation]
        return solution
    def Global_MA_LS(self, offspring): #끝까지 갔을 때 성능ㅇ 좋으면 변경
        coin = random.randint(1, 1)
        scheduling_seq ,factory_seq ,routing_seq, fittness, critical_machine, critical_operation = offspring
        if coin == 1:
            routing_seq2 = copy.deepcopy(routing_seq)
            machine = "M"+str(critical_machine)
            random.shuffle(critical_operation)
            factory = self.machine_factory(critical_machine)
            job = scheduling_seq[critical_operation[0]]
            k=0
            for i in range(critical_operation[0]+1):
                if job == scheduling_seq[i]:
                    k+=1
            job_op = "j"+str(job)+str(k)
            best_fittness = 10000
            best_routing = copy.deepcopy(routing_seq)
            for j in self.job_of_factory[factory]:
                if j != machine and self.p_table[j].loc[job_op] != 0 : 
                    routing_seq2[critical_operation[0]] = j
                    fittness2, critical_machine2,critical_operation2 = self.get_fittness(scheduling_seq, routing_seq2)
                    if fittness2 < best_fittness:
                        best_fittness = fittness2
                        best_routing = copy.deepcopy(routing_seq2)
            if best_fittness == 10000:
                solution = [scheduling_seq, factory_seq, routing_seq, fittness, critical_machine,critical_operation]
            else:
                fittness, critical_machine,critical_operation = self.get_fittness(scheduling_seq, best_routing)
                self.reset()
                solution = [scheduling_seq, factory_seq, routing_seq, fittness, critical_machine,critical_operation]
        return solution
    def Local_MA_LS(self, offspring): #local하게 좋은 
        coin = random.randint(1, 1)
        scheduling_seq ,factory_seq ,routing_seq, fittness, critical_machine = offspring
        if coin == 1:
            routing_seq2 = copy.deepcopy(routing_seq)
            critical_machine_index = []
            machine = "M"+str(critical_machine)#과부하 걸린 기계
            for i in range(40): #과부하 걸린 기계의 순서를 저장해놓음
                if machine == routing_seq[i]:
                    critical_machine_index.append(i) #[1,13,20,21,24,30,38]
            random.shuffle(critical_machine_index)
            factory = self.machine_factory(critical_machine)    
            job = scheduling_seq[critical_machine_index[0]] # 24의 7이었음
            k=0
            for i in range(critical_machine_index[0]+1): #24까지 돌면서 7이 몇번나오나
                if job == scheduling_seq[i]:
                    k+=1
            job_op = "j"+str(job)+str(k) #72
            for i in range(critical_machine_index[0]):
                operation = scheduling_seq[i]
                operation2 = "j" + str(operation) + str(self.job_preOperation[str(operation)])
                self.job_preOperation[str(operation)] += 1
                machine = routing_seq[i]
                time = max(self.machine_endTime[machine] ,self.job_endTime["j"+str(job)])
                df_sorted = self.p_table[machine]
                p_time = df_sorted.loc[operation2]
                df2_sorted = self.s_table[operation2[:2]]
                setup_time=df2_sorted.loc[self.machine_prejob[machine]]
                start = time
                end = start+setup_time+p_time
                self.machine_endTime[machine]=end #기계의 끝나는 시간 설정
                self.job_endTime[job]=end #job의 끝나는 시간 설정
                self.machine_prejob[machine] = "j"+str(job) #현재 어떤 machine에서 어떤 job을 수행했는지 기록
            machine = self.least_time_machine(job, job_op, factory_seq)
            routing_seq2[critical_machine_index[0]] = machine
            self.reset()
            fittness, critical_machine = self.get_fittness(scheduling_seq, routing_seq2)
            solution = [scheduling_seq, factory_seq, routing_seq2, fittness, critical_machine]
        return solution
    def replacement_operator2(self, population, offsprings):
        result_population = population[:]
        for i in range(5):
            fitness = offsprings[i][3]
            for j in range(100):
                if population[j][3] > fitness:
                    result_population[j] = copy.deepcopy(offsprings[i])
                    break
        return result_population
    def sort_population(self, population):
        population.sort(key=lambda x:x[3],reverse=False)
        # todo: fitness를 기준으로 population을 내림차순 정렬하고 반환
        return population
    def selection_operater(self, population):
        mom_ch = []
        dad_ch = []
        fitness_population=copy.deepcopy(population) #리스트 복사
        total_score = 0 #유전자의 적합도 총합 계산
        sum_fitness = 0 #0부터 적합도들의 합을 더함 이것이 k보다 클 시 그 해를 선택
        #유전자별 적합도를 새로운 리스트에 입력
        for i in range(len(population)):
            fitness_population[i][3]=abs((population[-1][3]-population[i][3])-(population[0][3]-population[-1][3])/(3-1))
        #(Cw-Ci)+(Cb-Cw)/(k-1)
        
        #유전자의 적합도 총합계산
        for i in range(len(fitness_population)):
            total_score=fitness_population[i][3]+total_score
        total_score=int(total_score)
        #랜덤함수 호출        
        k=random.uniform(0, total_score)
        #룰렛 돌리기
        for i in range(len(fitness_population)):
            sum_fitness=sum_fitness+fitness_population[i][3]
            if sum_fitness>=k:
                dad_ch=population[i]
                fitness_population.pop(i)
                break
        #적합도 총합을 다시계산해줌, 더했던 적합도들의 합도 다시 계산
        sum_fitness=0
        total_score=0
        for i in range(len(fitness_population)):
            total_score=fitness_population[i][3]+total_score
        total_score=int(total_score)
        #랜덤함수 실행
        k2=random.uniform(0, total_score)
        #룰렛 돌리기
        for i in range(len(fitness_population)):
            sum_fitness=sum_fitness+fitness_population[i][3]
            if sum_fitness>=k2:
                mom_ch=population[i]
                break
        return mom_ch, dad_ch
    def replacement_operator(self, population, offsprings):
        # todo: 생성된 자식해들(offsprings)을 이용하여 기존 해집단(population)의 해를 대치하여 새로운 해집단을 return
        """
        세대형 유전 알고리즘 사용
        해집단 내에서 가장 품질이 낮은 해를 대치하는 방법 사용(엘리티즘)
        """
        result_population = []
        for i in range(5):
            population.pop()
        for i in range(5):
            population.append(offsprings[i])
        result_population=population[:]
        return result_population
    def print_average_fitness(self, population):
        # todo: population의 평균 fitness를 출력
        population_average_fitness = 0
        total_population=0 
        for i in range(100):
            total_population=total_population+population[i][3]
        population_average_fitness=total_population/100
        print("population 평균 fitness: {}".format(population_average_fitness))
        return population_average_fitness
    def search(self):
        generation = 0  # 현재 세대 수
        population = [] # 해집단
        offsprings = [] # 자식해집단
        all_list = []                    
        for i in range(100):
            scheduling_seq = self.init_OS()
            routing_seq_f = self.init_FA()
            if i < 70:
                routing_seq_m = self.init_MA_Random(scheduling_seq, routing_seq_f)
            else:
                routing_seq_m = self.init_MA_LS(scheduling_seq, routing_seq_f)
            self.reset()
            fittness,critical_machine = self.get_fittness(scheduling_seq,routing_seq_m)
            self.reset()
            population.append([scheduling_seq,routing_seq_f,routing_seq_m,fittness,critical_machine])
        population = self.sort_population(population)
        result=[]
        while True:
            offsprings = []
            count_end=0 #동일 갯수
            for i in range(5):
                mom_ch, dad_ch = self.selection_operater(population)
                offspring = self.MOX_operator(dad_ch, mom_ch)
                self.reset()
                coin = random.random()
                if coin < self.eps2:
                    solution = self.Local_MA_LS(offspring)
                    self.reset()
                elif self.eps2 < coin < self.eps3:
                    solution = self.Two_Point_OS_LS(offspring)
                    self.reset()
                else:
                    coin2 = random.randint(0,1)
                    if coin2 == 0 :
                        solution = self.Local_FA_Move_LS(offspring)
                        self.reset()
                    else:
                        solution = self.Local_FA_Swap_LS(offspring)
                        self.reset()
                offsprings.append(solution)
            if generation > 100:
                population = self.replacement_operator2(population, offsprings)
            else:
                population = self.replacement_operator(population, offsprings)
            self.reset()
            population = self.sort_population(population)
            print('현재 세대',generation, '최고 해', population[0][3])
            avg = self.print_average_fitness(population)
            result.append([population[0][3],avg])
            self.anneal_eps()
            generation=generation+1
            pop_list=[]
            for i in range(100):
                pop_list.append(population[i][3])
            a = Counter(pop_list)
            end_number_list=[]
            for k in a.values():
                end_number_list.append(k)
            end_number = max(end_number_list)
            if end_number >= 90:
                self.reset()
                print(generation)
                break
            if generation == 1000:
                print(generation)
                break
        
        return population, result
total_makespan=0
makespan_list = []
starttime = datetime.now()
for i in range(20):
    if __name__ == "__main__":
        jaya = JAYA_FJSP(params)
        population,result = jaya.search()
        print(population[0])
        jaya.gannt_chart(population[0])
        total_makespan += population[0][3]
        makespan_list.append(population[0][3])
endtime = datetime.now()
elapsed_time = endtime-starttime
total_elapsed_time = 0
total_elapsed_time += elapsed_time.total_seconds()
print("100회 실행 평균 최적 makespan: ", total_makespan/20)
print("100회 평균 총 걸린시간: ", total_elapsed_time/20)
print("5회 최소값", min(makespan_list))
print("5회 최댓값", max(makespan_list))