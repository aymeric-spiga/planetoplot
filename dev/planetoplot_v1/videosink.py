import numpy as np
import subprocess
from os import system

class VideoSink(object) :

    def __init__( self, size, filename="output", rate=10, byteorder="bgra") :
            self.size = size
            cmdstring  = ('mencoder',
                    '/dev/stdin',
                    '-demuxer', 'rawvideo',
                    '-rawvideo', 'w=%i:h=%i'%size[::-1]+":fps=%i:format=%s"%(rate,byteorder),
                    '-o', filename+'_raw.avi',
                    '-nosound',
                    '-ovc', 'x264',
                    '-msglevel', 'all=-1'
#                    '-ovc', 'lavc',
                    )
            self.p = subprocess.Popen(cmdstring, stdin=subprocess.PIPE, shell=False)

    def run(self, image) :
            #assert image.shape[0:2] == self.size
            self.p.stdin.write(image.tostring())

    def first_pass(self,filename="output",quality=False,rate=10) :


             bitrate="7200"
             if quality:bitrate="50000"
             cmdstring  = ('mencoder',
                     'movie_raw.avi',
                     '-of', 'avi',
                     '-nosound',
                     '-ofps',str(rate),
                     '-ovc', 'x264',
                     '-x264encopts', 'subq=1:frameref=1:bitrate='+bitrate+':bframes=1:pass=1',
                     '-vf', 'scale=1280:720',
                     '-o', filename+'.avi'
                    )
             ## si second_pass '-o', filename+'_first.264'
             ## ... et seconde ligne: filename+'_raw.avi'
             subprocess.call(cmdstring,shell=False)

    def second_pass(self,filename="output",quality=False,rate=10) :
             bitrate="7200"
             if quality:bitrate="50000"
             cmdstring  = ('mencoder',
                     filename+'_first.264',
                     '-of', 'rawvideo',
                     '-nosound',
                     '-ofps',str(rate),
                     '-ovc', 'x264',
                     '-x264encopts', 'subq=6:frameref=5:bitrate='+bitrate+':me=umh:partitions=all:bframes=1:me_range=16:cabac:weightb:deblock:pass=2',
                     '-vf', 'scale=1280:720',
                     '-o', filename+'.avi'
                    )
             system('rm -f '+filename+'_raw.avi')
             subprocess.call(cmdstring,shell=False)
             system('rm -f '+filename+'_first.264')

    def close(self) :
            self.p.stdin.close()
            self.p.wait()

