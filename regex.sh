import re

with open("platemap_id_out.csv", "w") as outfile:
    with open("chemicals.txt", "w") as chem_file:
        failed = []
        regex = r"^(\d*),\"?(.*?)\s+((?:\d|\.){1,5}\s?(?:g/mL|mg/mL|mg|g|l|mL|ug))\"?(?:\t| )*(.*)?$"
        for line in open("platemap_id.csv", "r", encoding='cp1252').readlines():
            match = re.match(regex, line)
            if match == None:
                print(f"failed to match {line}")
                failed.append(line)
            else:
                groups = match.groups()
                three = groups[3].replace('\"','')
                newline = f"\"{groups[1]}\",{groups[2]},\"{three}\""
                outfile.write("success,"+newline + "\n")
                chem_file.write(groups[1].replace("\"","") + "\n")
        for line in failed:
            outfile.write("failed," + line)
            chem_file.write("failed on " + line)