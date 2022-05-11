'''
This python script is to generate buildin kernel source file,
if you change or add the kernel source *.cl, please run this script and
generate new cl_kernel_source.h
'''
import os
import shutil

def extension(file_path):
    return os.path.splitext(file_path)[1]

def main():
    files = os.listdir("./")
    files = [file for file in files if extension(file) == ".cl"]

    print("static const char *PROGRAM_SOURCE[] = {")
    currentLine = None
    for file in files:
        with open(file, "r") as src:
            while True:
                nextLine = src.readline()
                if currentLine is None:
                    currentLine = nextLine
                    nextLine = src.readline()

                if currentLine == "" or nextLine == "":
                    currentLine = "\"" + currentLine.rstrip() + "\\n\"" 
                    print(currentLine)
                    break
                else:
                    currentLine = "\"" + currentLine.rstrip() + "\\n\"," 
                    print(currentLine)
                    currentLine = nextLine;

                #if line == "":
                #    break
                #else: 
                #    line = "\"" + line.rstrip() + "\\n\"," 
                #    print(line)
    print("};")

if __name__ == "__main__":
    main()
