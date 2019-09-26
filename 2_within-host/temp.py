import re



s = ["001","120","078","101","008","123"]

for date in s:
    print(date)
    if re.search("^0+", date):
        date = re.sub(r"^0+(\d*)",r"\1", date)
        print(date)