# UCSC_Tracks
UCSC tracks source trackhub files

# Aging Mouse Brain
```python
upload_to_figshare(input_path="/ceph/gale-1/qzeng/AmbData/CellType_Allc/F.CellType.Age/*/*.CGN-both.frac.bw",
					   dataset_title='AgingMouseBrain', description='UCSC track hub',
					   token=None,output="figshare.tsv",rewrite=False,
					   mc_types=['CGN'],count_types=['frac'])

# prepare trackhub info
import os,sys
import pandas as pd
import glob

input_path="/ceph/gale-1/qzeng/AmbData/CellType_Allc/F.CellType.Age/*/*.CGN-both.frac.bw"
short_label='Aging Mouse Brain'
output="trackhub_info.txt"
    
input_path = os.path.abspath(os.path.expanduser(input_path))
output = os.path.abspath(os.path.expanduser(output))
if "*" not in input_path and os.path.isdir(input_path):
    input_files = [os.path.join(input_path, file) for file in os.listdir(input_path)]
else:
    input_files = glob.glob(input_path)
    
R=[]
for file_path in input_files:
    # print(input_file)
    file = os.path.basename(file_path)
    names = file.replace('.bw', '').split('.')
    ct, age, context, frac_cov = names
    mc_type=context.split('-')[0]
    R.append([ct,age,mc_type,frac_cov,file])
df=pd.DataFrame(R,columns=['parent','name','mc_type','count_type','file_name'])
df['ParentShortLabel']=df.parent.map(str)+'-'+df.mc_type.map(str)+'-'+df.count_type.map(str)
df['ParentLongLabel'] = short_label+ ' '+df.parent.map(str) + ' ' + df.mc_type.map(str)
df['TrackName']=df.ParentShortLabel+'_'+df.name
df['TrackShortLabel']=df.name
df['TrackLongLabel'] = df.ParentLongLabel.map(str) + ' '+df.TrackShortLabel
df['visibility'] = df.count_type.apply(lambda x:'dense' if x == 'frac' else 'hide')  # hide cov in default
df.to_csv(output,sep='\t',index=False)


# make trackhub
make_trackhub(
	trackhub_info="trackhub_info.txt",
	figshare_mapping="figshare.tsv",hub_name="AgingMouseBrain",
	short_label='Aging Mouse Brain',
    long_label='Aging Mouse Brain -Single Cell DNA Methylation',
	genome="mm10",email="qzeng@salk.edu",
	mc_types=['CGN','CHN'],count_types=['frac','cov'],
	max_mch=0.06)
```